#include "batchselect.h"

using namespace std;
using namespace seal;
using namespace seal::util;

void sample_poly_normal(
    shared_ptr<UniformRandomGenerator> prng, const EncryptionParameters &parms, uint64_t *destination,
    double noise_standard_deviation, double noise_max_deviation)
{
    auto &coeff_modulus = parms.coeff_modulus();
    size_t coeff_modulus_size = coeff_modulus.size();
    size_t coeff_count = parms.poly_modulus_degree();

    RandomToStandardAdapter engine(prng);
    ClippedNormalDistribution dist(
        0, noise_standard_deviation, noise_max_deviation);

    SEAL_ITERATE(iter(destination), coeff_count, [&](auto &I) {
        int64_t noise = static_cast<int64_t>(dist(engine));
        uint64_t flag = static_cast<uint64_t>(-static_cast<int64_t>(noise < 0));
        SEAL_ITERATE(
            iter(StrideIter<uint64_t *>(&I, coeff_count), coeff_modulus), coeff_modulus_size,
            [&](auto J) { *get<0>(J) = static_cast<uint64_t>(noise) + (flag & get<1>(J).value()); });
    });
}

void add_poly_error(
    size_t count,
    shared_ptr<UniformRandomGenerator> prng, const SEALContext::ContextData &context_data, uint64_t *destination,
    double noise_standard_deviation, double noise_max_deviation)
{
    auto &parms = context_data.parms();
    auto &coeff_modulus = parms.coeff_modulus();
    size_t coeff_modulus_size = coeff_modulus.size();
    size_t coeff_count = parms.poly_modulus_degree();

    Pointer<uint64_t> temp = allocate_poly(poly_modulus_degree, coeff_modulus_size, MemoryManager::GetPool());

    PolyIter destination_iter(destination, coeff_count, coeff_modulus_size);
    RNSIter temp_iter(temp.get(), coeff_count);

    for (size_t i = 0; i < count; ++i) {
        uint64_t *desti = destination+(coeff_modulus_size * coeff_count * i);
        sample_poly_normal(prng, parms, temp.get(), noise_standard_deviation, noise_max_deviation);
        ntt_negacyclic_harvey(temp_iter, coeff_modulus_size, context_data.small_ntt_tables());
        add_poly_coeffmod(destination_iter[i], temp_iter, coeff_modulus_size, coeff_modulus, destination_iter[i]);
    }
}

void sample_poly_uniform(
    shared_ptr<UniformRandomGenerator> prng, const EncryptionParameters &parms, uint64_t *destination)
{
    // Extract encryption parameters
    auto &coeff_modulus = parms.coeff_modulus();
    size_t coeff_modulus_size = coeff_modulus.size();
    size_t coeff_count = parms.poly_modulus_degree();
    size_t dest_byte_count = mul_safe(coeff_modulus_size, coeff_count, sizeof(uint64_t));

    constexpr uint64_t max_random = static_cast<uint64_t>(0xFFFFFFFFFFFFFFFFULL);

    // Fill the destination buffer with fresh randomness
    prng->generate(dest_byte_count, reinterpret_cast<seal_byte *>(destination));

    for (size_t j = 0; j < coeff_modulus_size; j++)
    {
        auto &modulus = coeff_modulus[j];
        uint64_t max_multiple = max_random - barrett_reduce_64(max_random, modulus) - 1;
        transform(destination, destination + coeff_count, destination, [&](uint64_t rand) {
            // This ensures uniform distribution
            while (rand >= max_multiple)
            {
                prng->generate(sizeof(uint64_t), reinterpret_cast<seal_byte *>(&rand));
            }
            return barrett_reduce_64(rand, modulus);
        });
        destination += coeff_count;
    }
}

// Destination needs to be a poly array of size len_n*len_m.
// The elements are ordered as a[1]*b[1],...,a[1]*b[len_b],a[2]*b[1],... (i.e., the matrix in a row-wise order)
void outer_product(PolyIter a, size_t len_a, PolyIter b, size_t len_b, PolyIter destination, const vector<Modulus> &coeff_modulus) {
    SEAL_ITERATE(iter(a, seq_iter(0)), len_a, [&](const tuple<RNSIter,uint64_t> &I) {
        SEAL_ITERATE(iter(b, destination + get<1>(I)*len_b), len_b, [&](const tuple<RNSIter,RNSIter> &J) {
             dyadic_product_coeffmod(get<0>(I), get<0>(J), coeff_modulus.size(), coeff_modulus, get<1>(J));
        });
    });
}

void vector_constant_product(PolyIter a, size_t len_a, RNSIter b, PolyIter destination, const vector<Modulus> &coeff_modulus) {
    auto deg = a.poly_modulus_degree();
    SEAL_ITERATE(iter(a, destination), len_a, [&](const tuple<RNSIter,RNSIter> &I) {
        dyadic_product_coeffmod(get<0>(I), b, coeff_modulus.size(), coeff_modulus, get<1>(I));
    });
}

void inner_product(PolyIter a, PolyIter b, size_t len, RNSIter destination, const vector<Modulus> &coeff_modulus) {
    auto temp(allocate_zero_poly(destination.poly_modulus_degree(), coeff_modulus.size(), MemoryPoolHandle::Global()));
    RNSIter tempIter(temp.get(), destination.poly_modulus_degree());
    set_zero_poly(destination.poly_modulus_degree(), coeff_modulus.size(), destination);
    SEAL_ITERATE(iter(a, b), len, [&](const tuple<RNSIter,RNSIter> &I) {
        dyadic_product_coeffmod(get<0>(I), get<1>(I), coeff_modulus.size(), coeff_modulus, tempIter);
        add_poly_coeffmod(tempIter, destination, coeff_modulus.size(), coeff_modulus, destination);
    });
}

void multiply_g(RNSIter a, PolyIter destination, const SEALContext::ContextData &context_data) {
    const EncryptionParameters &parms = context_data.parms();
    size_t poly_modulus_degree = parms.poly_modulus_degree();
    const vector<Modulus> &coeff_modulus = parms.coeff_modulus();
    size_t coeff_modulus_size = coeff_modulus.size();

    set_poly(a, poly_modulus_degree, coeff_modulus_size, destination);
    for (size_t i = 1; i < m; ++i) {
        multiply_poly_scalar_coeffmod(destination[i-1], coeff_modulus_size, g, coeff_modulus, destination[i]);
    }
}

void decompose_g(RNSIter y, PolyIter destination, const SEALContext::ContextData &context_data) {
    const EncryptionParameters &parms = context_data.parms();
    size_t poly_modulus_degree = parms.poly_modulus_degree();
    const vector<Modulus> &coeff_modulus = parms.coeff_modulus();
    size_t coeff_modulus_size = coeff_modulus.size();
    auto ntt_tables = context_data.small_ntt_tables();
    
    Pointer<uint64_t> y_composed = allocate_poly(poly_modulus_degree, coeff_modulus_size, MemoryManager::GetPool());
    RNSIter y_composed_iter(y_composed.get(), poly_modulus_degree);

    set_poly(*y, poly_modulus_degree, coeff_modulus_size, y_composed.get());
    inverse_ntt_negacyclic_harvey(y_composed_iter, coeff_modulus_size, ntt_tables); // inverse NTT
    context_data.rns_tool()->base_q()->compose_array(y_composed.get(), poly_modulus_degree, MemoryManager::GetPool()); // combine the two mod values into a single integers

    //cout << *y_composed.get() << " " << *(y_composed.get() + 1) << "\n";

    SEAL_ITERATE(destination, m, [&](const RNSIter &I) {
        // take mod g, and divide by g:
        //int counter = 0;
        SEAL_ITERATE(iter(StrideIter<uint64_t*>(y_composed, coeff_modulus_size), StrideIter<uint64_t*>((*I).ptr(), coeff_modulus_size)), poly_modulus_degree, [&](const tuple<uint64_t*,uint64_t*> &J) {
            //if (counter++ == 0) cout << *get<0>(J) << " " << *(get<0>(J)+1) << "\n";
            SEAL_DIVIDE_UINT128_UINT64(get<0>(J), g, get<1>(J));
            swap(*(get<0>(J)), *(get<1>(J)));
            swap(*(get<0>(J)+1), *(get<1>(J)+1));
        });
        context_data.rns_tool()->base_q()->decompose_array((*I).ptr(), poly_modulus_degree, MemoryManager::GetPool()); // back into mod form
        ntt_negacyclic_harvey(I, coeff_modulus_size, ntt_tables); // forward NTT
        //cout << counter++ << ":\n";
        //cout << poly_to_hex_string(I, poly_modulus_degree, 1) << "\n";//, pool);
        //cout << poly_to_hex_string(I + poly_modulus_degree, poly_modulus_degree, 1) << "\n";//, pool);
        //cout << poly_to_hex_string(y_decomposed.get() + (2*m+1)*poly_modulus_degree, poly_modulus_degree, 1) << "\n";//, pool);
    });
}




/**
Generates the public parameters.
*/
void LHE::setup() {
    const EncryptionParameters &parms = context_data_.parms();
    size_t poly_modulus_degree = parms.poly_modulus_degree();
    size_t coeff_modulus_size = parms.coeff_modulus().size();

    data_a_ = allocate_poly_array(w, poly_modulus_degree, coeff_modulus_size, MemoryManager::GetPool());
    SEAL_ITERATE(PolyIter(data_a_.get(), poly_modulus_degree, coeff_modulus_size), w, [&](const RNSIter &I) {
        sample_poly_uniform(prng, parms, I);
    });
    // no need to convert a into NTT, as we may assume that it was already sampled in NTT form
}

/**
Takes m1, and generates s1 and ct1 accordingly.
*/
void LHE::enc1(Pointer<uint64_t> &m1) {
    const EncryptionParameters &parms = context_data_.parms();
    size_t poly_modulus_degree = parms.poly_modulus_degree();
    const vector<Modulus> &coeff_modulus = parms.coeff_modulus();
    size_t coeff_modulus_size = coeff_modulus.size();

    data_s1_ = allocate_poly_array(m, poly_modulus_degree, coeff_modulus_size, MemoryManager::GetPool());
    data_ct1_ = allocate_poly_array(w*m, poly_modulus_degree, coeff_modulus_size, MemoryManager::GetPool());
    Pointer<uint64_t> temp = allocate_poly_array(m, poly_modulus_degree, coeff_modulus_size, MemoryManager::GetPool());

    PolyIter a_iter(data_a_.get(), poly_modulus_degree, coeff_modulus_size);
    PolyIter s1_iter(data_s1_.get(), poly_modulus_degree, coeff_modulus_size);
    PolyIter ct1_iter(data_ct1_.get(), poly_modulus_degree, coeff_modulus_size);
    PolyIter m1_iter(m1.get(), poly_modulus_degree, coeff_modulus_size);
    PolyIter temp_iter(temp.get(), poly_modulus_degree, coeff_modulus_size);

    //cout << "step 1" << endl;

    SEAL_ITERATE(s1_iter, m, [&](const RNSIter &I) {
        sample_poly_uniform(prng, parms, I);
    });
    // as for a, we just interprete s as polynomials in NTT form

    //cout << "step 2" << endl;

    outer_product(a_iter, w, s1_iter, m, ct1_iter, coeff_modulus);

    //cout << "step 3" << endl;

    for (size_t i = 0; i < w; ++i) {
        //cout << "step 5" << endl;
        multiply_g(m1_iter[i], temp_iter, context_data_);
        add_poly_coeffmod(ct1_iter + (i*m), temp_iter, m, coeff_modulus, ct1_iter + (i*m));
        //cout << "step 6" << endl;
    }
    //cout << "step 7" << endl;

    add_poly_error(w*m, prng, context_data_, data_ct1_.get(), noise_small_standard_deviation, noise_small_max_deviation);
}

/**
Takes m2, and generates s2 and ct2 accordingly.
*/
void LHE::enc2(Pointer<uint64_t> &m2) {
    const EncryptionParameters &parms = context_data_.parms();
    size_t poly_modulus_degree = parms.poly_modulus_degree();
    const vector<Modulus> &coeff_modulus = parms.coeff_modulus();
    size_t coeff_modulus_size = coeff_modulus.size();

    data_s2_ = allocate_poly(poly_modulus_degree, coeff_modulus_size, MemoryManager::GetPool());
    data_ct2_ = allocate_poly_array(w, poly_modulus_degree, coeff_modulus_size, MemoryManager::GetPool());

    PolyIter a_iter(data_a_.get(), poly_modulus_degree, coeff_modulus_size);
    RNSIter s2_iter(data_s2_.get(), poly_modulus_degree);
    PolyIter ct2_iter(data_ct2_.get(), poly_modulus_degree, coeff_modulus_size);
    PolyIter m2_iter(m2.get(), poly_modulus_degree, coeff_modulus_size);

    sample_poly_uniform(prng, parms, s2_iter);

    vector_constant_product(a_iter, w, s2_iter, ct2_iter, coeff_modulus);
    add_poly_coeffmod(ct2_iter, m2_iter, w, coeff_modulus, ct2_iter);

    add_poly_error(w, prng, context_data_, data_ct2_.get(), noise_large_standard_deviation, noise_large_max_deviation);
}

/**
Takes y, and computes sk_y from s1, s2, and y.
*/
void LHE::keygen(Pointer<uint64_t> &y) {
    const EncryptionParameters &parms = context_data_.parms();
    size_t poly_modulus_degree = parms.poly_modulus_degree();
    const vector<Modulus> &coeff_modulus = parms.coeff_modulus();
    size_t coeff_modulus_size = coeff_modulus.size();

    data_sk_ = allocate_poly(poly_modulus_degree, coeff_modulus_size, MemoryManager::GetPool());
    Pointer<uint64_t> y_decomposed = allocate_poly_array(m, poly_modulus_degree, coeff_modulus_size, MemoryManager::GetPool());
    Pointer<uint64_t> temp = allocate_poly(poly_modulus_degree, coeff_modulus_size, MemoryManager::GetPool());

    PolyIter s1_iter(data_s1_.get(), poly_modulus_degree, coeff_modulus_size);
    RNSIter sk_iter(data_sk_.get(), poly_modulus_degree);
    RNSIter y_iter(y.get(), poly_modulus_degree);
    PolyIter y_decomposed_iter(y_decomposed.get(), poly_modulus_degree, coeff_modulus_size);
    RNSIter temp_iter(temp.get(), poly_modulus_degree);

    decompose_g(y_iter, y_decomposed_iter, context_data_);
    /* for (int i = 0; i < m; ++i) {
        cout << i << ":\n";
        cout << poly_to_hex_string(y_decomposed.get() + (2*m)*poly_modulus_degree, poly_modulus_degree, 1) << "\n";//, pool);
        cout << poly_to_hex_string(y_decomposed.get() + (2*m+1)*poly_modulus_degree, poly_modulus_degree, 1) << "\n";//, pool);
    } */

    // sk <- s2
    set_poly(data_s2_.get(), poly_modulus_degree, coeff_modulus_size, data_sk_.get());
    SEAL_ITERATE(iter(s1_iter, y_decomposed_iter), m, [&](const tuple<RNSIter,RNSIter> &I) {
        // sk += s1[i]*y[i]
        dyadic_product_coeffmod(get<0>(I), get<1>(I), coeff_modulus_size, coeff_modulus, temp_iter);
        add_poly_coeffmod(sk_iter, temp_iter, coeff_modulus_size, coeff_modulus, sk_iter);
    });
}

/**
Takes y, and computes mres from ct1, ct2, and y.
*/
Pointer<uint64_t>& LHE::dec(Pointer<uint64_t> &y) {
    const EncryptionParameters &parms = context_data_.parms();
    size_t poly_modulus_degree = parms.poly_modulus_degree();
    const vector<Modulus> &coeff_modulus = parms.coeff_modulus();
    size_t coeff_modulus_size = coeff_modulus.size();

    data_mres_ = allocate_poly_array(w, poly_modulus_degree, coeff_modulus_size, MemoryManager::GetPool());
    Pointer<uint64_t> y_decomposed = allocate_poly_array(m, poly_modulus_degree, coeff_modulus_size, MemoryManager::GetPool());
    Pointer<uint64_t> temp = allocate_poly_array(w, poly_modulus_degree, coeff_modulus_size, MemoryManager::GetPool());

    PolyIter a_iter(data_a_.get(), poly_modulus_degree, coeff_modulus_size);
    RNSIter sk_iter(data_sk_.get(), poly_modulus_degree);
    PolyIter ct1_iter(data_ct1_.get(), poly_modulus_degree, coeff_modulus_size);
    PolyIter ct2_iter(data_ct2_.get(), poly_modulus_degree, coeff_modulus_size);
    PolyIter mres_iter(data_mres_.get(), poly_modulus_degree, coeff_modulus_size);
    RNSIter y_iter(y.get(), poly_modulus_degree);
    PolyIter y_decomposed_iter(y_decomposed.get(), poly_modulus_degree, coeff_modulus_size);
    PolyIter temp_iter(temp.get(), poly_modulus_degree, coeff_modulus_size);

    decompose_g(y_iter, y_decomposed_iter, context_data_);

    // mres <- ct1 * y
    for (size_t i = 0; i < w; ++i) {
        inner_product(ct1_iter + i*m, y_decomposed_iter, m, mres_iter[i], coeff_modulus);
    }
    // mres += ct2
    add_poly_coeffmod(mres_iter, ct2_iter, w, coeff_modulus, mres_iter);

    // mres -= a*sk
    vector_constant_product(a_iter, w, sk_iter, temp_iter, coeff_modulus);
    sub_poly_coeffmod(mres_iter, temp_iter, w, coeff_modulus, mres_iter);

    return data_mres_;
}





/**
Generates the public parameters.
*/
void Lenc::setup() {
    const EncryptionParameters &parms = context_data_.parms();
    size_t poly_modulus_degree = parms.poly_modulus_degree();
    size_t coeff_modulus_size = parms.coeff_modulus().size();

    data_b_ = allocate_poly_array(2*m, poly_modulus_degree, coeff_modulus_size, MemoryManager::GetPool());
    SEAL_ITERATE(PolyIter(data_b_.get(), poly_modulus_degree, coeff_modulus_size), 2*m, [&](const RNSIter &I) {
        sample_poly_uniform(prng, parms, I);
    });
    // no need to convert a into NTT, as we may assume that it was already sampled in NTT form
}

/**
Takes s, and generates r and ct accordingly.
*/
Pointer<uint64_t>& Lenc::enc(Pointer<uint64_t> &s) {
    const EncryptionParameters &parms = context_data_.parms();
    size_t poly_modulus_degree = parms.poly_modulus_degree();
    const vector<Modulus> &coeff_modulus = parms.coeff_modulus();
    size_t coeff_modulus_size = coeff_modulus.size();

    data_r_ = allocate_poly_array(l*w, poly_modulus_degree, coeff_modulus_size, MemoryManager::GetPool());
    data_ct_ = allocate_poly_array(l*w*2*m, poly_modulus_degree, coeff_modulus_size, MemoryManager::GetPool());
    Pointer<uint64_t> temp = allocate_poly_array(m, poly_modulus_degree, coeff_modulus_size, MemoryManager::GetPool());

    PolyIter b_iter(data_b_.get(), poly_modulus_degree, coeff_modulus_size);
    PolyIter r_iter(data_r_.get(), poly_modulus_degree, coeff_modulus_size);
    PolyIter s_iter(s.get(), poly_modulus_degree, coeff_modulus_size);
    PolyIter ct_iter(data_ct_.get(), poly_modulus_degree, coeff_modulus_size);
    PolyIter temp_iter(temp.get(), poly_modulus_degree, coeff_modulus_size);

    for (size_t i = 0; i < l; ++i) {
        PolyIter cti_iter = ct_iter + i*2*m*w;
        outer_product(r_iter + i*w, w, b_iter, 2*m, cti_iter, coeff_modulus);
    }
    cerr << "Lenc encryption first step done\n";
    chrono::nanoseconds g = chrono::steady_clock::now() - chrono::steady_clock::now();
    chrono::nanoseconds add = chrono::steady_clock::now() - chrono::steady_clock::now();
    for (size_t i = 0; i < l; ++i) {
        PolyIter cti_iter = ct_iter + i*2*m*w;
        for (size_t j = 0; j < w; ++j) {
            PolyIter ctij_iter = cti_iter + j*2*m;
            if (j & (1 << (l-i-1))) ctij_iter = ctij_iter + m;

            RNSIter ri_iter = i == l-1 ? s_iter[j] : r_iter[(i+1)*w + j];

            auto begin = chrono::steady_clock::now();
            multiply_g(ri_iter, temp_iter, context_data_);
            g += chrono::steady_clock::now() - begin;

            auto begin2 = chrono::steady_clock::now();
            add_poly_coeffmod(ctij_iter, temp_iter, m, coeff_modulus, ctij_iter);
            add += chrono::steady_clock::now() - begin2;
        }
    }

    cerr << "time for g: " << g.count() << "\n";
    cerr << "time for add: " << add.count() << "\n";

    auto begin = chrono::steady_clock::now();
    add_poly_error(l*w*2*m, prng, context_data_, data_ct_.get(), noise_small_standard_deviation, noise_small_max_deviation);
    cerr << "time for error: " << (chrono::steady_clock::now() - begin).count() << "\n";

    return data_r_;
}

/**
Takes a, and computes the digest y_epsilon.
*/
Pointer<uint64_t>& Lenc::digest(Pointer<uint64_t> &a) {
    const EncryptionParameters &parms = context_data_.parms();
    size_t poly_modulus_degree = parms.poly_modulus_degree();
    const vector<Modulus> &coeff_modulus = parms.coeff_modulus();
    size_t coeff_modulus_size = coeff_modulus.size();

    data_digest_ = allocate_poly(poly_modulus_degree, coeff_modulus_size, MemoryManager::GetPool());
    data_tree_ = allocate_poly_array((2*w-1)*m, poly_modulus_degree, coeff_modulus_size, MemoryManager::GetPool());
    Pointer<uint64_t> temp = allocate_poly(poly_modulus_degree, coeff_modulus_size, MemoryManager::GetPool());

    PolyIter a_iter(a.get(), poly_modulus_degree, coeff_modulus_size);
    PolyIter b_iter(data_b_.get(), poly_modulus_degree, coeff_modulus_size);
    PolyIter tree_iter(data_tree_.get(), poly_modulus_degree, coeff_modulus_size);
    RNSIter temp_iter(temp.get(), poly_modulus_degree);

    for (size_t i = 0; i < w; ++i) {
        decompose_g(a_iter[i], tree_iter + (w-1+i)*m, context_data_);
    }
    
    for (size_t i = w-1; i-- > 0;) {
        //cout << "going to multiply with ";
        //for (int k = 0; k < 2*m; ++k) if (poly_to_hex_string(tree_iter + (2*i+1)*m + k, poly_modulus_degree, 1) != "0") cout << "1 "; else cout << "0 ";
        //cout << "\n";
        inner_product(b_iter, tree_iter + (2*i+1)*m, 2*m, temp_iter, coeff_modulus);
        negate_poly_coeffmod(temp_iter, coeff_modulus_size, coeff_modulus, temp_iter);
        if (i) { // we do not need the decomposition of the root
            decompose_g(temp_iter, tree_iter + i*m, context_data_);
            //cout << poly_to_hex_string(temp_iter, poly_modulus_degree, 1) << "\n";//, pool);
            /* cout << i << ": ";
            for (size_t j = 0; j < 2*w-1; ++j) {
                bool one = false;
                for (int k = 0; k < m; ++k) if (poly_to_hex_string(tree_iter + j*m + k, poly_modulus_degree, 1) != "0") one = true;
                if (one) cout << "1 ";
                else cout << "0 ";
            }
            cout << "\n"; */
        } else { // instead, we will store the digest separately
            set_poly(temp.get(), poly_modulus_degree, coeff_modulus_size, data_digest_.get());
        }
    }

    return data_digest_;
}

/**
Takes a, and computes delta from ct and y.
digest(a) needs to be called first!
*/
Pointer<uint64_t>& Lenc::eval(Pointer<uint64_t> &a) {
    const EncryptionParameters &parms = context_data_.parms();
    size_t poly_modulus_degree = parms.poly_modulus_degree();
    const vector<Modulus> &coeff_modulus = parms.coeff_modulus();
    size_t coeff_modulus_size = coeff_modulus.size();

    data_delta_ = allocate_poly_array(w, poly_modulus_degree, coeff_modulus_size, MemoryManager::GetPool());
    Pointer<uint64_t> temp = allocate_poly(poly_modulus_degree, coeff_modulus_size, MemoryManager::GetPool());

    PolyIter delta_iter(data_delta_.get(), poly_modulus_degree, coeff_modulus_size);
    PolyIter ct_iter(data_ct_.get(), poly_modulus_degree, coeff_modulus_size);
    PolyIter tree_iter(data_tree_.get(), poly_modulus_degree, coeff_modulus_size);
    RNSIter temp_iter(temp.get(), poly_modulus_degree);

    SEAL_ITERATE(iter(seq_iter(0), delta_iter), w, [&](const tuple<uint64_t,RNSIter> &I) {
        inner_product(ct_iter + (2*m)*get<0>(I), tree_iter + m, 2*m, get<1>(I), coeff_modulus);
        for (size_t j = 1; j < l; ++j) {
            inner_product(ct_iter + (2*m*w)*j + (2*m)*get<0>(I), tree_iter + (((get<0>(I) >> (l-j)) + (1 << j) - 1)*2 + 1)*m, 2*m, temp_iter, coeff_modulus);
            add_poly_coeffmod(get<1>(I), temp_iter, coeff_modulus_size, coeff_modulus, get<1>(I));
        }
    });
    negate_poly_coeffmod(delta_iter, w, coeff_modulus, delta_iter);
    
    return data_delta_;
}




void BatchSelect::setup() {
    auto begin = chrono::steady_clock::now();
    cerr << "Setup...\n";
    lhe.setup();
    lenc.setup();
    cerr << "Setup done in " << (chrono::steady_clock::now() - begin).count() << ".\n";
}

void BatchSelect::enc1(Pointer<uint64_t> &l1) { // l1 should have size w * poly_modulus_degree
    const EncryptionParameters &parms = context_data_.parms();
    size_t poly_modulus_degree = parms.poly_modulus_degree();
    const vector<Modulus> &coeff_modulus = parms.coeff_modulus();
    size_t coeff_modulus_size = coeff_modulus.size();

    Pointer<uint64_t> temp = allocate_zero_poly_array(w, poly_modulus_degree, coeff_modulus_size, MemoryManager::GetPool());

    PolyIter temp_iter(temp.get(), poly_modulus_degree, coeff_modulus_size);

    for (size_t i = 0; i < w; ++i) {
        set_uint(l1.get() + i*poly_modulus_degree, poly_modulus_degree, temp_iter[i]);
        multiply_poly_scalar_coeffmod(temp_iter[i][0], poly_modulus_degree, coeff_modulus[1].value(), coeff_modulus[0], temp_iter[i][0]);
    }

    auto begin = chrono::steady_clock::now();
    cerr << "Lenc encryption...\n";
    Pointer<uint64_t> &r = lenc.enc(temp);
    cerr << "Lenc encryption done in " << (chrono::steady_clock::now() - begin).count() << ".\n";

    begin = chrono::steady_clock::now();
    cerr << "LHE encryption 1...\n";
    lhe.enc1(r);
    cerr << "LHE encryption 1 done in " << (chrono::steady_clock::now() - begin).count() << ".\n";
}

void BatchSelect::enc2(Pointer<uint64_t> &l2) {
    const EncryptionParameters &parms = context_data_.parms();
    size_t poly_modulus_degree = parms.poly_modulus_degree();
    const vector<Modulus> &coeff_modulus = parms.coeff_modulus();
    size_t coeff_modulus_size = coeff_modulus.size();

    Pointer<uint64_t> temp = allocate_zero_poly_array(w, poly_modulus_degree, coeff_modulus_size, MemoryManager::GetPool());

    PolyIter temp_iter(temp.get(), poly_modulus_degree, coeff_modulus_size);

    for (size_t i = 0; i < w; ++i) {
        set_uint(l2.get() + i*poly_modulus_degree, poly_modulus_degree, temp_iter[i]);
        multiply_poly_scalar_coeffmod(temp_iter[i][0], poly_modulus_degree, coeff_modulus[1].value(), coeff_modulus[0], temp_iter[i][0]);
    }

    add_poly_error(w, prng, context_data_, temp.get(), noise_large_standard_deviation, noise_large_max_deviation);

    auto begin = chrono::steady_clock::now();
    cerr << "LHE encryption 2...\n";
    lhe.enc2(temp);
    cerr << "LHE encryption 2 done in " << (chrono::steady_clock::now() - begin).count() << ".\n";
}

void BatchSelect::keygen(Pointer<uint64_t> &y) {
    const EncryptionParameters &parms = context_data_.parms();
    size_t poly_modulus_degree = parms.poly_modulus_degree();
    const vector<Modulus> &coeff_modulus = parms.coeff_modulus();
    size_t coeff_modulus_size = coeff_modulus.size();

    Pointer<uint64_t> temp = allocate_zero_poly_array(w, poly_modulus_degree, coeff_modulus_size, MemoryManager::GetPool());

    PolyIter temp_iter(temp.get(), poly_modulus_degree, coeff_modulus_size);

    for (size_t i = 0; i < w; ++i) {
        for (size_t j = 0; j < poly_modulus_degree; ++j) {
            if (y[i*poly_modulus_degree+j]) {
                temp_iter[i][0][j] = 1;
            }
        }
        set_uint(temp_iter[i][0], poly_modulus_degree, temp_iter[i][1]);
        inverse_ntt_negacyclic_harvey(temp_iter[i][1], *context_data_.plain_ntt_tables());
        ntt_negacyclic_harvey(temp_iter[i][1], context_data_.small_ntt_tables()[1]);
    }

    auto begin = chrono::steady_clock::now();
    cerr << "Computing Lenc digest (server)...\n";
    Pointer<uint64_t> &digest = lenc.digest(temp);
    cerr << "Computing Lenc digest (server) done in " << (chrono::steady_clock::now() - begin).count() << ".\n";

    begin = chrono::steady_clock::now();
    cerr << "LHE keygen...\n";
    lhe.keygen(digest);
    cerr << "LHE keygen done in " << (chrono::steady_clock::now() - begin).count() << ".\n";
}

void BatchSelect::dec(Pointer<uint64_t> &y, Pointer<uint64_t> &out) {
    const EncryptionParameters &parms = context_data_.parms();
    size_t poly_modulus_degree = parms.poly_modulus_degree();
    const vector<Modulus> &coeff_modulus = parms.coeff_modulus();
    size_t coeff_modulus_size = coeff_modulus.size();

    Pointer<uint64_t> temp = allocate_zero_poly_array(w, poly_modulus_degree, coeff_modulus_size, MemoryManager::GetPool());

    PolyIter temp_iter(temp.get(), poly_modulus_degree, coeff_modulus_size);

    for (size_t i = 0; i < w; ++i) {
        for (size_t j = 0; j < poly_modulus_degree; ++j) {
            if (y[i*poly_modulus_degree+j]) {
                temp_iter[i][0][j] = 1;
            }
        }
        set_uint(temp_iter[i][0], poly_modulus_degree, temp_iter[i][1]);
        inverse_ntt_negacyclic_harvey(temp_iter[i][1], *context_data_.plain_ntt_tables());
        ntt_negacyclic_harvey(temp_iter[i][1], context_data_.small_ntt_tables()[1]);
    }

    auto begin = chrono::steady_clock::now();
    cerr << "Computing Lenc digest (client)...\n";
    Pointer<uint64_t> &digest = lenc.digest(temp);
    cerr << "Computing Lenc digest (client) done in " << (chrono::steady_clock::now() - begin).count() << ".\n";

    begin = chrono::steady_clock::now();
    cerr << "LHE decryption...\n";
    Pointer<uint64_t> &res = lhe.dec(digest);
    cerr << "LHE decryption done in " << (chrono::steady_clock::now() - begin).count() << ".\n";

    begin = chrono::steady_clock::now();
    cerr << "Lenc evaluation...\n";
    Pointer<uint64_t> &delta = lenc.eval(temp);
    cerr << "Lenc evaluation done in " << (chrono::steady_clock::now() - begin).count() << ".\n";

    PolyIter res_iter(res.get(), poly_modulus_degree, coeff_modulus.size());
    PolyIter delta_iter(delta.get(), poly_modulus_degree, coeff_modulus.size());
    sub_poly_coeffmod(res_iter, delta_iter, w, coeff_modulus, res_iter);

    //cout << res_iter[0][0][0] << " " << res_iter[0][1][0] << "\n";
    //add_poly_error(1, prng, context_data_, res.get(), noise_small_standard_deviation, noise_small_max_deviation);
    //cout << res_iter[0][0][0] << " " << res_iter[0][1][0] << "\n";

    //cout << context_data_.rns_tool()->t().value() << " " << context_data_.rns_tool()->q_last_mod_t() << " " << context_data_.rns_tool()->m_sk().value() << " " << context_data_.rns_tool()->m_tilde().value() << " " << context_data_.rns_tool()->inv_q_last_mod_t() << " " << context_data_.rns_tool()->gamma().value() << "\n";
    //cout << context_data_.rns_tool()->base_q()->inv_punctured_prod_mod_base_array()->operand << " " << context_data_.rns_tool()->base_q()->inv_punctured_prod_mod_base_array()->quotient << "\n";
    MultiplyUIntModOperand inv = *context_data_.rns_tool()->base_q()->inv_punctured_prod_mod_base_array();
    for (size_t i = 0; i < w; ++i) {
        // Step 1: subtract the error from res_iter[i][0] (it is given in res_iter[i][1], but as a different modulus)
        inverse_ntt_negacyclic_harvey(res_iter[i][1], context_data_.small_ntt_tables()[1]);
        //if (i==0) cout << res_iter[0][0][0] << " " << res_iter[0][1][0] << "\n";
        // Step 1b: Now we need to convert to a different modulus (this is relevant whenever some coefficient is negative!)
        uint64_t modulus_value = coeff_modulus[1].value();
        uint64_t modulus_value_plaintext = coeff_modulus[0].value();
        SEAL_ITERATE(res_iter[i][1], poly_modulus_degree, [&](uint64_t &val) {
            //if (i==0) cout << "check " << val << " " << modulus_value/2 << "\n";
            if (val > modulus_value/2) val = modulus_value_plaintext - (modulus_value - val);
        });
        //if (i==0) cout << res_iter[0][0][0] << " " << res_iter[0][1][0] << "\n";
        ntt_negacyclic_harvey(res_iter[i][1], context_data_.small_ntt_tables()[0]);
        //if (i==0) cout << res_iter[0][0][0] << " " << res_iter[0][1][0] << "\n";
        sub_poly_coeffmod(res_iter[i][0], res_iter[i][1], poly_modulus_degree, coeff_modulus[0], res_iter[i][0]);
        // Step 2: remove the factor of Delta from res_iter[i][0] by multiplying with its inverse
        multiply_poly_scalar_coeffmod(res_iter[i][0], poly_modulus_degree, inv, coeff_modulus[0], res_iter[i][0]);
        // Step 3: copy output
        set_uint(res_iter[i], poly_modulus_degree, out.get() + i*poly_modulus_degree);
    }
    //cout << res_iter[0][0][0] << " " << res_iter[0][1][0] << "\n";
}