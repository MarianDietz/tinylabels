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

    BatchSelect bs(context_data, prng);
    bs.setup();

    FILE *f_pp = fopen("pp.bin", "wb");
    bs.save_pp(f_pp);
    fclose(f_pp);

    cout << "# Add's = " << counter_poly_add << "\n";
    cout << "# Sub's = " << counter_poly_sub << "\n";
    cout << "# Mult's = " << counter_poly_mult << "\n";
    cout << "# Scalar mult's = " << counter_poly_mult_scalar << "\n";
    cout << "# Neg's = " << counter_poly_negate << "\n";
    cout << "# Forward NTT's = " << counter_ntt_forward << "\n";
    cout << "# Inverse NTT's = " << counter_ntt_inverse << "\n";

    // Pointer<uint64_t> l1(allocate_zero_uint(w*poly_modulus_degree, MemoryManager::GetPool()));
    // Pointer<uint64_t> l2(allocate_zero_uint(w*poly_modulus_degree, MemoryManager::GetPool()));
    // Pointer<uint64_t> out(allocate_zero_uint(w*poly_modulus_degree, MemoryManager::GetPool()));
    // vector<bool> y(w*poly_modulus_degree, false);

    // for (size_t i = 0; i < w*poly_modulus_degree; ++i) {
    //     l1[i] = rand()%coeff_modulus[0].value();
    //     l2[i] = rand()%coeff_modulus[0].value();
    //     y[i] = rand()%2;
    // }

    // cout << coeff_modulus[0].value() << " " << coeff_modulus[1].value() << "\n";
    // bs.enc1(l1);
    // bs.enc2(l2);
    // bs.keygen(y);
    // bs.dec(y, out);

    // cerr << "Batch-select done.\n";

    // for (size_t i = 0; i < w*poly_modulus_degree; ++i) {
    //     uint64_t expected = (l1[i]*y[i] + l2[i]) % coeff_modulus[0].value();
    //     if (expected != out[i]) cout << "incorrect: " << i << "\n";
    // }

    // cerr << "Diffcheck done.\n";

    return 0;
}
