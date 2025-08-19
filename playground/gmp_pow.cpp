// montmul_compare_montgomery.cpp
#include "gmpxx.h"
#include <iostream>
#include <vector>
#include <chrono>
#include <cstdlib>
#include <cstring>

int main(int argc, char* argv[]) {
    if (argc != 2) {
        std::cerr << "Usage: " << argv[0] << " <number_of_integers>\n";
        return 1;
    }
    int ncount = std::atoi(argv[1]);
    if (ncount <= 0) {
        std::cerr << "Number of integers must be positive.\n";
        return 1;
    }

    const unsigned int BITS = 1024;

    gmp_randclass rng(gmp_randinit_default);
    rng.seed(static_cast<unsigned long>(time(nullptr)));

    // Generate random numbers (plain integers)
    std::vector<mpz_class> nums;
    nums.reserve(ncount);
    for (int i = 0; i < ncount; ++i) nums.emplace_back(rng.get_z_bits(BITS));

    // Generate random odd modulus
    mpz_class modulus = rng.get_z_bits(BITS);
    modulus |= mpz_class(1) << (BITS - 1);
    modulus |= 1;

    // -------------------------------
    // 1) High-level GMP timing: (a*b) % modulus
    // -------------------------------
    auto start_high = std::chrono::high_resolution_clock::now();
    for (int i = 0; i < ncount; ++i) {
        for (int j = 0; j < ncount; ++j) {
            mpz_class prod = (nums[i] * nums[j]) % modulus;
            (void)prod;
        }
    }
    auto end_high = std::chrono::high_resolution_clock::now();
    double total_ms_high = std::chrono::duration<double,std::milli>(end_high - start_high).count();

    long long num_mults = (long long)ncount * (long long)ncount;

    // gmp_printf ("%Zd\n", nums[0].get_mpz_t()); // Print first num
    std::cout << "=== High-level GMP (a*b % p) ===\n";
    std::cout << "Number of group multiplications: " << num_mults << "\n";
    std::cout << "Total time (loop only): " << total_ms_high << " ms\n";
    std::cout << "Average per multiplication: " << (total_ms_high / num_mults * 1000.0) << " µs\n\n";

    // -------------------------------
    // 2) Powers timing: (a^b) % modulus
    // -------------------------------
    //

    // Generate random 256 bit numbers (plain integers)
    std::vector<mpz_class> nums1;
    nums1.reserve(ncount);
    int n_bits_exponent = 160;
    // int n_bits_exponent = 256;
    for (int i = 0; i < ncount; ++i) nums1.emplace_back(rng.get_z_bits(n_bits_exponent));

    auto start_pow = std::chrono::high_resolution_clock::now();
    for (int i = 0; i < ncount; ++i) {
        for (int j = 0; j < ncount; ++j) {
            mpz_class pow_res = mpz_class(nums[i].get_mpz_t());
            mpz_powm(pow_res.get_mpz_t(), pow_res.get_mpz_t(), nums1[j].get_mpz_t(), modulus.get_mpz_t());
            (void)pow_res;
        }
    }
    auto end_pow = std::chrono::high_resolution_clock::now();
    auto total_ms_pow = std::chrono::duration<double,std::milli>(end_pow - start_pow).count();
    // Timing
    std::cout << "=== High-level GMP (a^b % p) ===\n";
    std::cout << "Total time (loop only): " << total_ms_pow << " ms\n";
    std::cout << "Average per power: " << (total_ms_pow / (ncount * ncount) * 1000.0) << " µs\n\n";

    return 0;
}
