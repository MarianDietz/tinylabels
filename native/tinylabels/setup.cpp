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

    auto begin = chrono::steady_clock::now();

    BatchSelect bs(context_data, prng);
    bs.setup();

    cout << "===================\n";
    cout << "Total time: " << time_str(chrono::steady_clock::now() - begin) << ".\n";

    print_statistics();

    FILE *f_pp = fopen("pp.bin", "wb");
    bs.save_pp(f_pp);
    fclose(f_pp);

    return 0;
}
