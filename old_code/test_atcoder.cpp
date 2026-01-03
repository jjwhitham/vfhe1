#include "convolution.hpp"
#include "gmpxx.h"

int main() {
    // mpz_class mod = 11;
    std::vector<mpz_class> a = {0, -1};
    std::vector<mpz_class> b = {5, 5};
    // atcoder::static_modint::set_m(mpz_class(17)); // Set the modulus for static_modint
    atcoder::static_modint::set_m(mpz_class(18014398509506561)); // Set the modulus for static_modint
    std::vector<mpz_class> conv = atcoder::convolution(a, b);
    std::cout << "Convolution result: ";
    for (const auto& val : conv)
        gmp_printf("%Zd ", val.get_mpz_t());
    std::cout << "\n";
    return 0;
}