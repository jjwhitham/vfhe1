// montmul_compare_montgomery.cpp
#include <gmpxx.h>
#include <iostream>
#include <vector>
#include <chrono>
#include <cstdlib>
#include <cstring>

typedef __uint128_t u128;
// Compute modular inverse of m0 (least significant limb of modulus) mod 2^GMP_NUMB_BITS
static mp_limb_t inv_mod_base(mp_limb_t m0) {
    mp_limb_t inv = 1;
    for (int i = 0; i < 6; ++i) { // Newton iteration (enough for 64-bit limbs)
        u128 t = (u128)m0 * inv;
        mp_limb_t prod = (mp_limb_t)t;
        mp_limb_t two_minus = (mp_limb_t)(2 - prod);
        u128 next = (u128)inv * two_minus;
        inv = (mp_limb_t)next;
    }
    return (mp_limb_t)(0 - inv);
}

// Montgomery multiply: res = (a * b * R^{-1}) mod m
static void montgomery_mul(mp_limb_t* res,
                           const mp_limb_t* a,
                           const mp_limb_t* b,
                           const mp_limb_t* m,
                           mp_size_t n,
                           mp_limb_t m0inv_neg) {
    std::vector<mp_limb_t> t(2*n + 1);
    std::fill(t.begin(), t.end(), 0);

    mpn_mul(t.data(), a, n, b, n); // t = a * b

    for (mp_size_t i = 0; i < n; ++i) {
        mp_limb_t u = (mp_limb_t)((u128)t[i] * (u128)m0inv_neg);
        mp_limb_t carry1 = mpn_addmul_1(t.data() + i, m, n, u);
        mp_limb_t c = mpn_add_1(t.data() + i + n, t.data() + i + n, 1, carry1);
        if (c) { // propagate carry
            mp_size_t pos = i + n + 1;
            while (c && pos < static_cast<long int>(t.size())) {
                u128 sum = (u128)t[pos] + c;
                t[pos] = (mp_limb_t)sum;
                c = (mp_limb_t)(sum >> (sizeof(mp_limb_t)*8));
                ++pos;
            }
        }
    }

    // Copy out t[n .. n+n-1] into res
    std::memcpy(res, t.data() + n, n * sizeof(mp_limb_t));

    // If res >= m, subtract once
    if (mpn_cmp(res, m, n) >= 0) {
        mpn_sub_n(res, res, m, n);
    }
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

    const unsigned int BITS = 3000;
    const mp_size_t n_limbs = (BITS + GMP_NUMB_BITS - 1) / GMP_NUMB_BITS;

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

    // Export modulus to limbs
    std::vector<mp_limb_t> m_limbs(n_limbs);
    size_t count;
    mpz_export(m_limbs.data(), &count, -1, sizeof(mp_limb_t), 0, 0, modulus.get_mpz_t());

    // Precompute Montgomery constant (m0inv_neg)
    mp_limb_t m0inv_neg = inv_mod_base(m_limbs[0]);

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

    // -------------------------------
    // Prepare Montgomery-domain representations
    // -------------------------------
    // R = base^{n_limbs} where base = 2^{GMP_NUMB_BITS}
    mpz_class R = mpz_class(1) << (n_limbs * GMP_NUMB_BITS);
    mpz_class R_mod = R % modulus; // R mod p

    // Convert each num to Montgomery domain: a_mont = (a * R) mod p
    std::vector<mpz_class> nums_mont;
    nums_mont.reserve(ncount);

    auto start_conv = std::chrono::high_resolution_clock::now();
    for (int i = 0; i < ncount; ++i) {
        nums_mont.emplace_back((nums[i] * R_mod) % modulus);
    }
    auto end_conv = std::chrono::high_resolution_clock::now();
    double total_ms_conv = std::chrono::duration<double,std::milli>(end_conv - start_conv).count();

    // Export Montgomery numbers to limb arrays
    std::vector<std::vector<mp_limb_t>> nums_limbs_mont;
    nums_limbs_mont.resize(ncount, std::vector<mp_limb_t>(n_limbs, 0));
    for (int i = 0; i < ncount; ++i) {
        mpz_export(nums_limbs_mont[i].data(), &count, -1, sizeof(mp_limb_t), 0, 0, nums_mont[i].get_mpz_t());
    }

    // Also export a limb-array for "1" in Montgomery domain to convert back (montone = 1)
    std::vector<mp_limb_t> one_limbs(n_limbs, 0);
    mpz_class one = 1;
    mpz_class one_mont = (one * R_mod) % modulus; // representation of 1 in Montgomery domain
    mpz_export(one_limbs.data(), &count, -1, sizeof(mp_limb_t), 0, 0, one_mont.get_mpz_t());

    std::vector<mp_limb_t> res_limbs(n_limbs);

    // -------------------------------
    // 2) Montgomery timing (in-domain multiplications)
    // -------------------------------
    auto start_mont = std::chrono::high_resolution_clock::now();
    for (int i = 0; i < ncount; ++i) {
        for (int j = 0; j < ncount; ++j) {
            // inputs are in Montgomery domain; montgomery_mul preserves Montgomery domain
            montgomery_mul(res_limbs.data(),
                           nums_limbs_mont[i].data(),
                           nums_limbs_mont[j].data(),
                           m_limbs.data(),
                           n_limbs,
                           m0inv_neg);
        }
    }
    auto end_mont = std::chrono::high_resolution_clock::now();
    double total_ms_mont = std::chrono::duration<double,std::milli>(end_mont - start_mont).count();

    // Optional: convert one result back to regular representation (res * 1 * R^{-1} => normal)
    // This is done by montgomery_mul(res, one_limbs) once (example only)
    montgomery_mul(res_limbs.data(), res_limbs.data(), one_limbs.data(), m_limbs.data(), n_limbs, m0inv_neg);
    // If you want to interpret res_limbs as mpz, you'd mpz_import it.

    long long num_mults = (long long)ncount * (long long)ncount;

    std::cout << "=== High-level GMP (a*b % p) ===\n";
    std::cout << "Total time (loop only): " << total_ms_high << " ms\n";
    std::cout << "Average per multiplication: " << (total_ms_high / num_mults * 1000.0) << " µs\n\n";

    std::cout << "=== Montgomery approach ===\n";
    std::cout << "Conversion to Montgomery (n numbers): " << total_ms_conv << " ms\n";
    std::cout << "Montgomery loop total (in-domain multiplies): " << total_ms_mont << " ms\n";
    std::cout << "Total Mont (conversion + loop): " << (total_ms_conv + total_ms_mont) << " ms\n";
    std::cout << "Average per multiplication (Mont loop only): " << (total_ms_mont / num_mults * 1000.0) << " µs\n";
    std::cout << "Average per multiplication (including conversion): " << ((total_ms_conv + total_ms_mont) / num_mults * 1000.0) << " µs\n\n";

    if (total_ms_mont > 0.0) {
        std::cout << "Speedup (high-level / mont loop only): " << (total_ms_high / total_ms_mont) << "x\n";
        std::cout << "Speedup (high-level / mont total incl conv): " << (total_ms_high / (total_ms_conv + total_ms_mont)) << "x\n";
    }

    return 0;
}
