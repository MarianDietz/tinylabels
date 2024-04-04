#include "mytests.h"

using namespace std;
using namespace seal;
using namespace seal::util;

size_t w = 10;
size_t m = 1;

constexpr double noise_standard_deviation = 3.2;
constexpr double noise_max_deviation = 128 * noise_standard_deviation;

void sample_poly_normal(
    shared_ptr<UniformRandomGenerator> prng, const EncryptionParameters &parms, uint64_t *destination)
{
    auto coeff_modulus = parms.coeff_modulus();
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

void sample_poly_uniform(
    shared_ptr<UniformRandomGenerator> prng, const EncryptionParameters &parms, uint64_t *destination)
{
    // Extract encryption parameters
    auto coeff_modulus = parms.coeff_modulus();
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
// The elements are ordered as a[1]*b[1],...,a[1]*b[len_m],a[2]*b[1],... (i.e., the matrix in a row-wise order)
void outer_product(PolyIter a, size_t len_a, PolyIter b, size_t len_b, PolyIter destination, vector<Modulus> &coeff_modulus) {
    SEAL_ITERATE(iter(a, seq_iter(0)), len_a, [&](auto I) {
        SEAL_ITERATE(iter(b, destination + get<1>(I)*len_b), len_b, [&](auto J) {
             dyadic_product_coeffmod(get<0>(I), get<0>(J), coeff_modulus.size(), coeff_modulus, get<1>(J));
        });
    });
}

void vector_constant_product(PolyIter a, size_t len_a, MultiplyUIntModOperand constant, PolyIter destination, vector<Modulus> &coeff_modulus) {
    auto deg = a.poly_modulus_degree();
    SEAL_ITERATE(iter(a, destination), len_a, [&](auto I) {
        SEAL_ITERATE(iter(get<0>(I), coeff_modulus, get<1>(I)), coeff_modulus.size(), [&](auto J) {
            multiply_poly_scalar_coeffmod(get<0>(J), deg, constant, get<1>(J), get<2>(J));
        });
    });
}

void inner_product(PolyIter a, PolyIter b, size_t len, RNSIter destination, vector<Modulus> &coeff_modulus) {
    auto temp(allocate_poly(destination.poly_modulus_degree(), coeff_modulus.size(), MemoryPoolHandle::Global()));
    RNSIter tempIter(temp.get(), destination.poly_modulus_degree());
    SEAL_ITERATE(iter(a, b), len, [&](auto I) {
        dyadic_product_coeffmod(get<0>(I), get<1>(I), coeff_modulus.size(), coeff_modulus, tempIter);
        add_poly_coeffmod(tempIter, destination, coeff_modulus.size(), coeff_modulus, destination);
    });
    // TODO may have to deallocate temp unless SEAL does this automatically
}

int main()
{
    cout << "Microsoft SEAL version: " << SEAL_VERSION << endl;

    EncryptionParameters parms(scheme_type::bfv);
    size_t poly_modulus_degree = 8192;
    parms.set_poly_modulus_degree(poly_modulus_degree);
    parms.set_plain_modulus(PlainModulus::Batching(poly_modulus_degree, 43));

    vector<Modulus> coeff_modulus = CoeffModulus::Create(poly_modulus_degree, {44});
    coeff_modulus.push_back(parms.plain_modulus());
    parms.set_coeff_modulus(coeff_modulus);

    SEALContext context(parms);
    print_parameters(context);
    cout << endl;

    auto &context_data = *context.get_context_data(parms.parms_id());
    auto ntt_tables = context_data.small_ntt_tables();

    cout << "first coeff modulus: " << parms.coeff_modulus()[0].value() << endl;

    // auto qualifiers = context.first_context_data()->qualifiers();
    // cout << qualifiers.parameter_error_message() << endl;

    BatchEncoder batch_encoder(context);
    size_t slot_count = batch_encoder.slot_count();
    vector<uint64_t> pod_matrix(slot_count, 0ULL);
    pod_matrix[0] = 0;

    Plaintext plain_matrix;
    batch_encoder.encode(pod_matrix, plain_matrix);

    cout << "is ntt: " << plain_matrix.is_ntt_form() << endl;
    cout << "data: " << *((uint64_t*) plain_matrix.data()) << endl;

    auto prng = UniformRandomGeneratorFactory::DefaultFactory()->create();
    MemoryPoolHandle pool = MemoryManager::GetPool(mm_prof_opt::mm_force_new, true);

    // auto noise(allocate_poly(poly_modulus_degree, coeff_modulus.size(), pool));
    // sample_poly_normal(prng, parms, noise.get());
    // for (int i = 0; i < coeff_modulus.size(); ++i) {
    //     ntt_negacyclic_harvey(noise.get() + i * poly_modulus_degree, ntt_tables[i]);
    // }

    // RNSIter noiseIter(noise.get(), poly_modulus_degree);
    // modulo_poly_coeffs(noiseIter, coeff_modulus.size(), iter(coeff_modulus), noiseIter);

    auto m1(allocate_zero_poly_array(w, poly_modulus_degree, coeff_modulus.size(), pool));
    PolyIter m1Iter(m1.get(), poly_modulus_degree, coeff_modulus.size());
    auto m2(allocate_zero_poly_array(w, poly_modulus_degree, coeff_modulus.size(), pool));
    PolyIter m2Iter(m2.get(), poly_modulus_degree, coeff_modulus.size());
    auto y(allocate_zero_poly(poly_modulus_degree, coeff_modulus.size(), pool));
    RNSIter yIter(y.get(), poly_modulus_degree);

    m1[3] = 2;
    m2[5] = 1;
    y[0] = 3;

    auto a(allocate_poly_array(w, poly_modulus_degree, coeff_modulus.size(), pool));
    PolyIter aIter(a.get(), poly_modulus_degree, coeff_modulus.size());
    SEAL_ITERATE(aIter, w, [&](auto &I) {
        sample_poly_uniform(prng, parms, I);
    });
    // no need to convert a into NTT, as we may assume that it was already sampled in NTT form

    // s[0]..s[m-1] denote s_1, s[m] denotes s_2
    auto s(allocate_poly_array(m + 1, poly_modulus_degree, coeff_modulus.size(), pool));
    PolyIter sIter(s.get(), poly_modulus_degree, coeff_modulus.size());
    SEAL_ITERATE(sIter, w, [&](auto &I) {
        sample_poly_uniform(prng, parms, I);
    });
    // as for a, we just interprete s as polynomials in NTT form

    auto ct1(allocate_poly_array(w*m, poly_modulus_degree, coeff_modulus.size(), pool));
    PolyIter ct1Iter(ct1.get(), poly_modulus_degree, coeff_modulus.size());
    outer_product(aIter, w, sIter, m, ct1Iter, coeff_modulus);
    add_poly_coeffmod(ct1Iter, m1Iter, w, coeff_modulus, ct1Iter);

    auto ct2(allocate_poly_array(w, poly_modulus_degree, coeff_modulus.size(), pool));
    PolyIter ct2Iter(ct2.get(), poly_modulus_degree, coeff_modulus.size());
    outer_product(aIter, w, sIter+m, 1, ct2Iter, coeff_modulus);
    add_poly_coeffmod(ct2Iter, m2Iter, w, coeff_modulus, ct2Iter);

    auto ctres(allocate_poly_array(w, poly_modulus_degree, coeff_modulus.size(), pool));
    PolyIter ctresIter(ctres.get(), poly_modulus_degree, coeff_modulus.size());
    outer_product(ct1Iter, w, PolyIter(y.get(), poly_modulus_degree, coeff_modulus.size()), 1, ctresIter, coeff_modulus);
    add_poly_coeffmod(ctresIter, ct2Iter, w, coeff_modulus, ctresIter);

    auto sk(allocate_poly(poly_modulus_degree, coeff_modulus.size(), pool));
    RNSIter skIter(sk.get(), poly_modulus_degree);
    inner_product(sIter, PolyIter(y.get(), poly_modulus_degree, coeff_modulus.size()), 1, skIter, coeff_modulus);
    add_poly_coeffmod(skIter, *(sIter+m), coeff_modulus.size(), coeff_modulus, skIter);

    auto aSk(allocate_poly_array(w, poly_modulus_degree, coeff_modulus.size(), pool));
    PolyIter aSkIter(aSk.get(), poly_modulus_degree, coeff_modulus.size());
    outer_product(aIter, w, PolyIter(sk.get(), poly_modulus_degree, coeff_modulus.size()), 1, aSkIter, coeff_modulus);

    auto mres(allocate_poly_array(w, poly_modulus_degree, coeff_modulus.size(), pool));
    PolyIter mresIter(mres.get(), poly_modulus_degree, coeff_modulus.size());
    sub_poly_coeffmod(ctresIter, aSkIter, w, coeff_modulus, mresIter);

    cout << poly_to_hex_string(mres.get(), poly_modulus_degree, 1) << "\n";//, pool);
    cout << poly_to_hex_string(mres.get() + poly_modulus_degree, poly_modulus_degree, 1) << "\n";//, pool);

    return 0;
}
