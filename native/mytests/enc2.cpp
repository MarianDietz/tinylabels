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

    cout << "Plaintext modulus: " << coeff_modulus[0].value() << "\n";
    cout << "===================\n";

    BatchSelect bs(context_data, prng);

    FILE *f_pp = fopen("pp.bin", "rb");
    bs.read_pp(f_pp);
    fclose(f_pp);

    Pointer<uint64_t> l2(allocate_zero_uint(w*poly_modulus_degree, MemoryManager::GetPool()));
    ifstream f_l2;
    f_l2.open("l2.txt");
    for (size_t i = 0; i < w*poly_modulus_degree; ++i) {
        f_l2 >> l2[i];
    }
    f_l2.close();

    auto begin = chrono::steady_clock::now();

    bs.enc2(l2);

    cout << "===================\n";
    cout << "Total time: " << time_str(chrono::steady_clock::now() - begin) << ".\n";
    print_statistics();

    FILE *f_st2 = fopen("st2.bin", "wb");
    bs.save_st2(f_st2);
    fclose(f_st2);

    FILE *f_ct2 = fopen("ct2.bin", "wb");
    bs.save_ct2(f_ct2);
    fclose(f_ct2);

    return 0;
}
