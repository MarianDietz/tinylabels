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

    cout << "Plaintext modulus: " << coeff_modulus[0].value() << "\n";
    cout << "===================\n";

    auto prng = UniformRandomGeneratorFactory::DefaultFactory()->create();

    BatchSelect bs(context_data, prng);

    FILE *f_pp = fopen("pp.bin", "rb");
    bs.read_pp(f_pp);
    fclose(f_pp);

    FILE *f_st1 = fopen("st1.bin", "rb");
    bs.read_st1(f_st1);
    fclose(f_st1);

    FILE *f_st2 = fopen("st2.bin", "rb");
    bs.read_st2(f_st2);
    fclose(f_st2);

    Pointer<uint64_t> y(allocate_zero_uint(w*poly_modulus_degree, MemoryManager::GetPool()));
    ifstream f_y;
    f_y.open("y.txt");
    for (size_t i = 0; i < w*poly_modulus_degree; ++i) {
        f_y >> y[i];
    }
    f_y.close();

    auto begin = chrono::steady_clock::now();

    bs.keygen(y);
    
    cout << "===================\n";
    cout << "Total time: " << time_str(chrono::steady_clock::now() - begin) << ".\n";
    print_statistics();

    FILE *f_sk = fopen("sk.bin", "wb");
    bs.save_sk(f_sk);
    fclose(f_sk);

    return 0;
}
