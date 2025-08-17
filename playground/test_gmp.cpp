#include <gmpxx.h>
#include <iostream>
#include <vector>
#include <chrono>
#include <cstdlib>
#include <cstring>

// Compute modular inverse of m0 (least significant limb of modulus) mod 2^GMP_NUMB_BITS
static mp_limb_t inv_mod_base(mp_limb_t m0) {
    mp_limb_t inv = 1;
    for (int i = 0; i < 6; ++i) { // Newton iteration
        unsigned __int128 t = (unsigned __int128)m0 * inv;
        mp_limb_t prod = (mp_limb_t)t;
        mp_limb_t two_minus = (mp_limb_t)(2 - prod);
        unsigned __int128 next = (unsigned __int128)inv * two_minus;
        inv = (mp_limb_t)next;
    }
    return (mp_limb_t)(0 - inv); // negate for Montgomery usage
}

// Montgomery multiply: res = (a*b*R^{-1}) mod m
static void montgomery_mul(mp_limb_t* res,
                           const mp_limb_t* a,
                           const mp_limb_t* b,
                           const mp_limb_t* m,
                           mp_size_t n,
                           mp_limb_t m0inv_neg) {
    std::vector<mp_limb_t> t(2*n + 1, 0);
    mpn_mul(t.data(), a, n, b, n); // t = a*b

    for (mp_size_t i = 0; i < n; ++i) {
        mp_limb_t u = t[i] * m0inv_neg; // truncated mod base
        mp_limb_t carry1 = mpn_addmul_1(t.data() + i, m, n, u);
        mp_limb_t c = mpn_add_1(t.data() + i + n, t.data() + i + n, 1, carry1);
        if (c) { // propagate carry
            mp_size_t pos = i + n + 1;
            while (c && pos < t.size()) {
                unsigned __int128 sum = (unsigned __int128)t[pos] + c;
                t[pos] = (mp_limb_t)sum;
                c = (mp_limb_t)(sum >> (sizeof(mp_limb_t)*8));
                ++pos;
            }
        }
    }

    // Copy out t[n..n+n-1] to res
    std::memcpy(res, t.data() + n, n * sizeof(mp_limb_t));

    // Conditional subtraction
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

    // Generate random numbers
    std::vector<mpz_class> nums;
    nums.reserve(ncount);
    for (int i = 0; i < ncount; ++i) {
        nums.emplace_back(rng.get_z_bits(BITS));
    }

    // Generate random odd modulus
    mpz_class modulus = rng.get_z_bits(BITS);
    modulus |= mpz_class(1) << (BITS - 1);
    modulus |= 1;

    // Export modulus to limbs
    std::vector<mp_limb_t> m_limbs(n_limbs);
    size_t count;
    mpz_export(m_limbs.data(), &count, -1, sizeof(mp_limb_t), 0, 0, modulus.get_mpz_t());

    // Montgomery precompute
    mp_limb_t m0inv_neg = inv_mod_base(m_limbs[0]);

    // Export all nums to limb arrays
    std::vector<std::vector<mp_limb_t>> nums_limbs;
    nums_limbs.resize(ncount, std::vector<mp_limb_t>(n_limbs, 0));
    for (int i = 0; i < ncount; ++i) {
        mpz_export(nums_limbs[i].data(), &count, -1, sizeof(mp_limb_t), 0, 0, nums[i].get_mpz_t());
    }

    std::vector<mp_limb_t> res_limbs(n_limbs);

    // -------------------------------
    // High-level GMP timing
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
    // Montgomery timing
    // -------------------------------
    auto start_mont = std::chrono::high_resolution_clock::now();
    for (int i = 0; i < ncount; ++i) {
        for (int j = 0; j < ncount; ++j) {
            montgomery_mul(res_limbs.data(),
                           nums_limbs[i].data(),
                           nums_limbs[j].data(),
                           m_limbs.data(),
                           n_limbs,
                           m0inv_neg);
        }
    }
    auto end_mont = std::chrono::high_resolution_clock::now();
    double total_ms_mont = std::chrono::duration<double,std::milli>(end_mont - start_mont).count();

    long long num_mults = (long long)ncount * (long long)ncount;

    std::cout << "=== High-level GMP (a*b % m) ===\n";
    std::cout << "Total time: " << total_ms_high << " ms\n";
    std::cout << "Average: " << (total_ms_high / num_mults * 1000.0) << " µs\n";

    std::cout << "\n=== Low-level Montgomery ===\n";
    std::cout << "Total time: " << total_ms_mont << " ms\n";
    std::cout << "Average: " << (total_ms_mont / num_mults * 1000.0) << " µs\n";

    std::cout << "\nSpeedup factor: " << (total_ms_high / total_ms_mont) << "x\n";

    return 0;
}
