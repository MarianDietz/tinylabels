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

    FILE *f_ct1 = fopen("ct1.bin", "rb");
    bs.read_ct1(f_ct1);
    fclose(f_ct1);

    FILE *f_ct2 = fopen("ct2.bin", "rb");
    bs.read_ct2(f_ct2);
    fclose(f_ct2);

    FILE *f_sk = fopen("sk.bin", "rb");
    bs.read_sk(f_sk);
    fclose(f_sk);

    Pointer<uint64_t> y(allocate_zero_uint(w*poly_modulus_degree, MemoryManager::GetPool()));

    FILE *f_y = fopen("y.bin", "rb");
    fread(y.get(), 8, w*poly_modulus_degree, f_y);
    fclose(f_y);

    Pointer<uint64_t> out(allocate_zero_uint(w*poly_modulus_degree, MemoryManager::GetPool()));

    bs.dec(y, out);

    FILE *f_output = fopen("output.bin", "wb");
    fwrite(out.get(), 8, w*poly_modulus_degree, f_output);
    fclose(f_output);

    cout << "# Add's = " << counter_poly_add << "\n";
    cout << "# Sub's = " << counter_poly_sub << "\n";
    cout << "# Mult's = " << counter_poly_mult << "\n";
    cout << "# Scalar mult's = " << counter_poly_mult_scalar << "\n";
    cout << "# Neg's = " << counter_poly_negate << "\n";
    cout << "# Forward NTT's = " << counter_ntt_forward << "\n";
    cout << "# Inverse NTT's = " << counter_ntt_inverse << "\n";

    return 0;
}
