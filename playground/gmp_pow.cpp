// montmul_compare_montgomery.cpp
#include "gmpxx.h"
#include <iostream>
#include <vector>
#include <chrono>
#include <cstdlib>
#include <cstring>

mpz_class mod_(mpz_class val, mpz_class q) {
    val %= q;
    if (val < 0) {
        val += q;
    }
    return val;
}

// mpz_class pow_(mpz_class base, mpz_class power, mpz_class FIELD_MODULUS, mpz_class GROUP_MODULUS) {
mpz_class pow_(mpz_class base, mpz_class power, mpz_class GROUP_MODULUS) {
    // power = mod_(power, FIELD_MODULUS);
    // base = mod_(base, GROUP_MODULUS);
    mpz_class result = 1;
    while (power > 0) {
        bool is_power_odd = (power % 2) == 1;
        if (is_power_odd)
            result = (result * base) % GROUP_MODULUS;
        power >>= 1;
        base = (base * base) % GROUP_MODULUS;
    }
    return result;
}

void test_mult_and_pow(int ncount, int n_bits_base, int n_bits_exponent) {
    gmp_randclass rng(gmp_randinit_default);
    rng.seed(static_cast<unsigned long>(time(nullptr)));

    // Generate random numbers (plain integers)
    std::vector<mpz_class> nums;
    nums.reserve(ncount);
    for (int i = 0; i < ncount; ++i) nums.emplace_back(rng.get_z_bits(n_bits_base));

    // Generate random odd modulus
    mpz_class modulus = rng.get_z_bits(n_bits_base);
    modulus |= mpz_class(1) << (n_bits_base - 1);
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
    std::cout << "ops/sec: " << (num_mults / total_ms_high * 1000.0) << "\n\n";

    // -------------------------------
    // 2) Powers timing: (a^b) % modulus
    // -------------------------------
    //

    // Generate random n_bit_ numbers (plain integers)
    std::vector<mpz_class> nums1;
    nums1.reserve(ncount);
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
    std::cout << "Average per power: " << (total_ms_pow / num_mults * 1000.0) << " µs\n\n";
    std::cout << "ops/sec: " << (num_mults / total_ms_pow * 1000.0) << "\n\n";

//     // -------------------------------
//     // 3) My pow timing: (a*b) % modulus
//     // -------------------------------
//     auto start_pow1 = std::chrono::high_resolution_clock::now();
//     for (int i = 0; i < ncount; ++i) {
//         for (int j = 0; j < ncount; ++j) {
//             mpz_class exp = pow_(nums[i], nums[j], modulus);
//             (void)exp;
//         }
//     }
//     auto end_pow1 = std::chrono::high_resolution_clock::now();
//     auto total_ms_pow1 = std::chrono::duration<double,std::milli>(end_pow1 - start_pow1).count();

//     std::cout << "=== My pow (a^b % p) ===\n";
//     std::cout << "Total time (loop only): " << total_ms_pow1 << " ms\n";
//     std::cout << "Average per power: " << (total_ms_pow1 / num_mults * 1000.0) << " µs\n\n";
//     std::cout << "ops/sec: " << (num_mults / total_ms_pow1 / 1000.0) << " µs\n\n";
}

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

    int n_bits_base = 2048;
    int n_bits_exponent = 224;
    test_mult_and_pow(ncount, n_bits_base, n_bits_exponent);
    return 0;
}
