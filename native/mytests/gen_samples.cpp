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
    cout << "Generating l1, l2, y, and expected (l1*y+l2)...\n";

    Pointer<uint64_t> l1(allocate_zero_uint(w*poly_modulus_degree, MemoryManager::GetPool()));
    Pointer<uint64_t> l2(allocate_zero_uint(w*poly_modulus_degree, MemoryManager::GetPool()));
    Pointer<uint64_t> out(allocate_zero_uint(w*poly_modulus_degree, MemoryManager::GetPool()));
    Pointer<uint64_t> y(allocate_zero_uint(w*poly_modulus_degree, MemoryManager::GetPool()));

    std::random_device rd;
    std::mt19937_64 e2(rd());
    std::uniform_int_distribution<uint64_t> dist(0, coeff_modulus[0].value()-1);
    std::uniform_int_distribution<uint64_t> bin(0, 1);

    for (size_t i = 0; i < w*poly_modulus_degree; ++i) {
        l1[i] = dist(e2);
        l2[i] = dist(e2);
        y[i] = bin(e2);
        out[i] = (l1[i]*y[i] + l2[i]) % coeff_modulus[0].value();
    }

    cout << "Done.\n";

    ofstream f_l1;
    f_l1.open("l1.txt");
    for (size_t i = 0; i < w*poly_modulus_degree; ++i) {
        f_l1 << l1[i] << "\n";
    }
    f_l1.close();

    ofstream f_l2;
    f_l2.open("l2.txt");
    for (size_t i = 0; i < w*poly_modulus_degree; ++i) {
        f_l2 << l2[i] << "\n";
    }
    f_l2.close();

    ofstream f_y;
    f_y.open("y.txt");
    for (size_t i = 0; i < w*poly_modulus_degree; ++i) {
        f_y << y[i] << "\n";
    }
    f_y.close();

    ofstream f_expected;
    f_expected.open("expected.txt");
    for (size_t i = 0; i < w*poly_modulus_degree; ++i) {
        f_expected << out[i] << "\n";
    }
    f_expected.close();

    return 0;
}
