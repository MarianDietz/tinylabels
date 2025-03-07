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
    ifstream f_y;
    f_y.open("y.txt");
    for (size_t i = 0; i < w*poly_modulus_degree; ++i) {
        f_y >> y[i];
    }
    f_y.close();

    auto begin = chrono::steady_clock::now();

    Pointer<uint64_t> out(allocate_zero_uint(w*poly_modulus_degree, MemoryManager::GetPool()));
    bs.dec(y, out);

    cout << "===================\n";
    cout << "Total time: " << time_str(chrono::steady_clock::now() - begin) << ".\n";
    print_statistics();

    ofstream f_out;
    f_out.open("output.txt");
    for (size_t i = 0; i < w*poly_modulus_degree; ++i) {
        f_out << out[i] << "\n";
    }
    f_out.close();

    return 0;
}
