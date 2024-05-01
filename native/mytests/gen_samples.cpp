#include "batchselect.h"

using namespace std;
using namespace seal;
using namespace seal::util;

int main()
{
    EncryptionParameters parms(scheme_type::onoff);
    parms.set_poly_modulus_degree(poly_modulus_degree);
    parms.set_plain_modulus(PlainModulus::Batching(poly_modulus_degree, 43));

    vector<Modulus> coeff_modulus = CoeffModulus::Create(poly_modulus_degree, {44});
    coeff_modulus.insert(coeff_modulus.begin(), parms.plain_modulus());
    parms.set_coeff_modulus(coeff_modulus);

    SEALContext context(parms);
    auto &context_data = *context.get_context_data(parms.parms_id());

    auto prng = UniformRandomGeneratorFactory::DefaultFactory()->create();

    Pointer<uint64_t> l1(allocate_zero_uint(w*poly_modulus_degree, MemoryManager::GetPool()));
    Pointer<uint64_t> l2(allocate_zero_uint(w*poly_modulus_degree, MemoryManager::GetPool()));
    Pointer<uint64_t> out(allocate_zero_uint(w*poly_modulus_degree, MemoryManager::GetPool()));
    Pointer<uint64_t> y(allocate_zero_uint(w*poly_modulus_degree, MemoryManager::GetPool()));

    for (size_t i = 0; i < w*poly_modulus_degree; ++i) {
        l1[i] = rand()%coeff_modulus[0].value();
        l2[i] = rand()%coeff_modulus[0].value();
        y[i] = rand()%2;
        out[i] = (l1[i]*y[i] + l2[i]) % coeff_modulus[0].value();
    }

    cout << coeff_modulus[0].value() << "\n";
    for (size_t i = 0; i < 5; ++i) cout << l1[i] << "*" << y[i] << " + " << l2[i] << " = " << out[i] << "\n";

    FILE *f_l1 = fopen("l1.bin", "wb");
    fwrite(l1.get(), 8, w*poly_modulus_degree, f_l1);
    fclose(f_l1);

    FILE *f_l2 = fopen("l2.bin", "wb");
    fwrite(l2.get(), 8, w*poly_modulus_degree, f_l2);
    fclose(f_l2);

    FILE *f_y = fopen("y.bin", "wb");
    fwrite(y.get(), 8, w*poly_modulus_degree, f_y);
    fclose(f_y);

    FILE *f_expected = fopen("expected.bin", "wb");
    fwrite(out.get(), 8, w*poly_modulus_degree, f_expected);
    fclose(f_expected);

    return 0;
}
