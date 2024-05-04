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

    SEALContext context(parms);
    auto &context_data = *context.get_context_data(parms.parms_id());

    auto prng = UniformRandomGeneratorFactory::DefaultFactory()->create();

    BatchSelect bs(context_data, prng);

    FILE *f_pp = fopen("pp.bin", "rb");
    bs.read_pp(f_pp);
    fclose(f_pp);

    Pointer<uint64_t> l1(allocate_zero_uint(w*poly_modulus_degree, MemoryManager::GetPool()));

    FILE *f_l1 = fopen("l1.bin", "rb");
    fread(l1.get(), 8, w*poly_modulus_degree, f_l1);
    fclose(f_l1);

    bs.enc1(l1);

    FILE *f_st1 = fopen("st1.bin", "wb");
    bs.save_st1(f_st1);
    fclose(f_st1);

    FILE *f_ct1 = fopen("ct1.bin", "wb");
    bs.save_ct1(f_ct1);
    fclose(f_ct1);

    cout << "# Add's = " << counter_poly_add << " (" << time_poly_add.count() << " ns)\n";
    cout << "# Sub's = " << counter_poly_sub << " (" << time_poly_sub.count() << " ns)\n";
    cout << "# Mult's = " << counter_poly_mult << " (" << time_poly_mult.count() << " ns)\n";
    cout << "# Scalar mult's = " << counter_poly_mult_scalar << " (" << time_poly_mult_scalar.count() << " ns)\n";
    cout << "# Forward NTT's = " << counter_ntt_forward << " (" << time_ntt_forward.count() << " ns)\n";
    cout << "# Inverse NTT's = " << counter_ntt_inverse << " (" << time_ntt_inverse.count() << " ns)\n";
    cout << "# Compose = " << counter_poly_compose << " (" << time_poly_compose.count() << " ns)\n";
    cout << "# Decompose = " << counter_poly_decompose << " (" << time_poly_decompose.count() << " ns)\n";

    return 0;
}
