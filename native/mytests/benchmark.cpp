#include "batchselect.h"

using namespace std;
using namespace seal;
using namespace seal::util;

int main()
{
    EncryptionParameters parms(scheme_type::onoff);
    parms.set_poly_modulus_degree(poly_modulus_degree);
    parms.set_plain_modulus(PlainModulus::Batching(poly_modulus_degree, mod_plaintext));

    vector<Modulus> coeff_modulus = CoeffModulus::Create(poly_modulus_degree, {mod_noise});
    coeff_modulus.insert(coeff_modulus.begin(), parms.plain_modulus());
    parms.set_coeff_modulus(coeff_modulus);

    cout << "moduli: " << coeff_modulus[0].value() << " " << coeff_modulus[1].value() << "\n";

    SEALContext context(parms);
    auto &context_data = *context.get_context_data(parms.parms_id());

    auto prng = UniformRandomGeneratorFactory::DefaultFactory()->create();
    
    Pointer<uint64_t> a = allocate_zero_poly_array(5000, poly_modulus_degree, coeff_modulus.size(), MemoryManager::GetPool());
    Pointer<uint64_t> b = allocate_zero_poly_array(5000, poly_modulus_degree, coeff_modulus.size(), MemoryManager::GetPool());
    Pointer<uint64_t> c = allocate_zero_poly_array(5000, poly_modulus_degree, coeff_modulus.size(), MemoryManager::GetPool());

    for (size_t i = 0; i < 5000; ++i) {
        sample_poly_uniform(prng, parms, a.get());
        sample_poly_uniform(prng, parms, b.get());
    }

    PolyIter a_iter(a.get(), poly_modulus_degree, coeff_modulus.size());
    PolyIter b_iter(b.get(), poly_modulus_degree, coeff_modulus.size());
    PolyIter c_iter(c.get(), poly_modulus_degree, coeff_modulus.size());

    auto begin = chrono::steady_clock::now();
    for (size_t i = 0; i < 5000; ++i) {
        add_poly_coeffmod(a_iter[i], b_iter[i], coeff_modulus.size(), coeff_modulus, c_iter[i]);
    }
    cerr << "10000 add's done in " << time_str(chrono::steady_clock::now() - begin) << ".\n";

    begin = chrono::steady_clock::now();
    for (size_t i = 0; i < 5000; ++i) {
        sub_poly_coeffmod(a_iter[i], b_iter[i], coeff_modulus.size(), coeff_modulus, c_iter[i]);
    }
    cerr << "10000 sub's done in " << time_str(chrono::steady_clock::now() - begin) << ".\n";

    begin = chrono::steady_clock::now();
    for (size_t i = 0; i < 5000; ++i) {
        dyadic_product_coeffmod(a_iter[i], b_iter[i], coeff_modulus.size(), coeff_modulus, c_iter[i]);
    }
    cerr << "10000 mult's done in " <<time_str(chrono::steady_clock::now() - begin) << ".\n";

    begin = chrono::steady_clock::now();
    for (size_t i = 0; i < 5000; ++i) {
        multiply_poly_scalar_coeffmod(a_iter[i], coeff_modulus.size(), i, coeff_modulus, c_iter[i]);
    }
    cerr << "10000 scalar mult's done in " << time_str(chrono::steady_clock::now() - begin) << ".\n";

    begin = chrono::steady_clock::now();
    for (size_t i = 0; i < 5000; ++i) {
        negate_poly_coeffmod(a_iter[i], coeff_modulus.size(), coeff_modulus, c_iter[i]);
    }
    cerr << "10000 neg's done in " << time_str(chrono::steady_clock::now() - begin) << ".\n";

    begin = chrono::steady_clock::now();
    for (size_t i = 0; i < 5000; ++i) {
        ntt_negacyclic_harvey(a_iter[i], coeff_modulus.size(), context_data.small_ntt_tables());
    }
    cerr << "10000 forward NTT's done in " << time_str(chrono::steady_clock::now() - begin) << ".\n";

    begin = chrono::steady_clock::now();
    for (size_t i = 0; i < 5000; ++i) {
        inverse_ntt_negacyclic_harvey(a_iter[i], coeff_modulus.size(), context_data.small_ntt_tables());
    }
    cerr << "10000 inverse NTT's done in " << time_str(chrono::steady_clock::now() - begin) << ".\n";

    print_statistics();

    return 0;
}
