#pragma once

#include <iostream>
#include <vector>
#include <cassert>
#include <cmath>
// #include <gmpxx.h>
#include <string>
#include <array>
#include <span>
#include "vfhe.h"
#include "shared.h"

using u128 = __uint128_t;
using u32 = unsigned int;
using vec_u128 = std::vector<u128>;

constexpr u32 POLY_SIZE = N_;
constexpr u128 FIELD_MOD = FIELD_MODULUS;
// constexpr u128 NTH_ROU = 2;
// constexpr u128 TWO_ROU = 3;
constexpr u128 GROUP_MOD = GROUP_MODULUS;
// p = k * q + 1
// constexpr u128 K = 6;
// constexpr u128 GENERATOR = 64;

using arr2n_u128 = std::array<u128, 2 * POLY_SIZE>;
using arr_u128 = std::array<u128, POLY_SIZE>;
using span_u128 = std::span<u128>;

u128 mod_sub(u128 x, u128 y) {
    return (x < y) ? ((x + FIELD_MOD) - y) : (x - y);
}

u128 pow_(u128 base, u128 power, u128 mod) {
    u128 result = 1;
    while (power > 0) {
        bool power_is_odd = (power & 1);
        if (power_is_odd)
            result = (result * base) % mod;
        base = (base * base) % mod;
        power >>= 1;
    }
    return result;
}

constexpr u128 pow_constexpr(u128 base, u128 power, u128 mod) {
    u128 result = 1;
    while (power > 0) {
        bool power_is_odd = (power & 1);
        if (power_is_odd)
            result = (result * base) % mod;
        base = (base * base) % mod;
        power >>= 1;
    }
    return result;
}

/* START ATCODER */
/* Code reproduced from atcoder (CC0 1.0 licence)
https://github.com/atcoder/ac-library/tree/master?tab=readme-ov-file */

// @param m `1 <= m`
// @return x mod m
constexpr long long safe_mod(long long x, long long m) {
    x %= m;
    if (x < 0) x += m;
    return x;
}

// @param n `0 <= n`
// @param m `1 <= m`
// @return `(x ** n) % m`
constexpr long long pow_mod_constexpr(long long x, long long n, int m) {
    if (m == 1) return 0;
    unsigned int _m = (unsigned int)(m);
    unsigned long long r = 1;
    unsigned long long y = safe_mod(x, m);
    while (n) {
        if (n & 1) r = (r * y) % _m;
        y = (y * y) % _m;
        n >>= 1;
    }
    return r;
}
// Reference:
// M. Forisek and J. Jancina,
// Fast Primality Testing for Integers That Fit into a Machine Word
// @param n `0 <= n`
constexpr bool is_prime_constexpr(int n) {
    if (n <= 1) return false;
    if (n == 2 || n == 7 || n == 61) return true;
    if (n % 2 == 0) return false;
    long long d = n - 1;
    while (d % 2 == 0) d /= 2;
    constexpr long long bases[3] = {2, 7, 61};
    for (long long a : bases) {
        long long t = d;
        long long y = pow_mod_constexpr(a, t, n);
        while (t != n - 1 && y != 1 && y != n - 1) {
            y = y * y % n;
            t <<= 1;
        }
        if (y != n - 1 && t % 2 == 0) {
            return false;
        }
    }
    return true;
}
/* END ATCODER */

constexpr bool are_q_and_N_legal(u128 q, u32 N) {
    bool both_legal = true;
    // DONE check q is prime
    both_legal = both_legal && is_prime_constexpr(q);
    // q % n == 1
    both_legal = both_legal && (q % N == 1);
    // q % 2n == 1
    both_legal = both_legal && (q % (2 * N) == 1);
    // log2(n) % 1 < 0.0001
    // or, pow = round(log2(n)); 2^pow == n;
    bool is_pow_2 = (N & (N - 1)) == 0;
    both_legal = both_legal && is_pow_2;
    return both_legal;
}

constexpr bool is_w_legal(u128 w, u32 n=POLY_SIZE) {
    // w is an nth ROU: w^n == 1
    bool w_pow_n_is_one = pow_constexpr(w, n, FIELD_MOD) == 1;
    if (!w_pow_n_is_one)
        return false;
    // w is a primitive ROU: w^i != 1, for all i = [1, n - 1]
    for (size_t i = 1; i < n; i++) {
        if (pow_constexpr(w, i, FIELD_MOD) == 1)
            return false;
    }
    return true;
}

constexpr bool is_psi_legal(u128 psi) {
    // psi is a 2n-th ROU: psi^2n == 1
    // psi is a primitive 2n-th ROU: psi^i != 1, for all i = [1, 2n - 1]
    u32 TWO_N = 2 * POLY_SIZE;
    return is_w_legal(psi, TWO_N);
}

constexpr arr_u128 get_rou_pows(u128 rou) {
    arr_u128 w_pows{};
    w_pows[0] = 1;
    for (size_t i = 1; i < POLY_SIZE; i++) {
        w_pows[i] = (w_pows[i - 1] * rou) % FIELD_MOD;
    }
    return w_pows;
}

constexpr bool is_w_pows_legal(const arr_u128& w_pows, u32 order=POLY_SIZE) {
    if (w_pows[0] != 1) return false;
    for (size_t i = 0; i < POLY_SIZE; i++) {
        if (pow_constexpr(w_pows[i], order, FIELD_MOD) != 1)
            return false;
    }
    return true;
}

constexpr bool is_psi_pows_legal(const arr_u128& psi_pows) {
    return is_w_pows_legal(psi_pows, 2 * POLY_SIZE);
}

// sorts the array in bit-reversed order
void bit_reverse(auto& a) {
    size_t n = a.size();
    size_t j = 0;
    for (size_t i = 1; i < n; i++) {
        size_t bit = n >> 1;
        while (j & bit) {
            j ^= bit;
            bit >>= 1;
        }
        j |= bit;
        if (i < j) {
            std::swap(a[i], a[j]);
        }
    }
}

void print_arr(auto& x) {
    std::cout << "x: ";
    for (auto& el : x)
        std::cout << static_cast<unsigned long int>(el) << " ";
    std::cout << "\n";
}

void ntt_recursive_(span_u128 x, u128 nth_root_unity) {
    size_t n = x.size();
    if (n == 1)
        return;
    size_t half = n / 2;
    vec_u128 w(half);
    w[0] = 1;
    for (size_t i = 1; i < half; i++)
        w[i] = (w[i - 1] * nth_root_unity) % FIELD_MOD;
    auto x_even = x.subspan(0, half);
    auto x_odd = x.subspan(half, half);
    nth_root_unity = pow_(nth_root_unity, 2, FIELD_MOD);
    ntt_recursive_(x_even, nth_root_unity);
    ntt_recursive_(x_odd, nth_root_unity);
    for (size_t i = 0; i < half; i++) {
        auto even = x_even[i];
        auto w_odd = (w[i] * x_odd[i]) % FIELD_MOD;
        x[i] = (even + w_odd) % FIELD_MOD;
        x[i + half] = mod_sub(even, w_odd);
    }
}

void ntt_iter_(auto& x, const arr_u128& w) {
    // const size_t n = x.size();
    constexpr size_t n = 2 * POLY_SIZE;
    // TODO alg that does no->bo fwd and bo->no rev could remove this step
    bit_reverse(x);

    for (size_t half = 1; half < n; half *= 2) {
        const size_t len_chunk = 2 * half;
        const size_t n_chunks = n / len_chunk;
        for (size_t chunk = 0; chunk < n_chunks; chunk++) {
            const size_t j = len_chunk * chunk;
            for (size_t i = 0; i < half; i++) {
                const auto even = x.get(i + j);
                const auto w_odd = (w[i * n_chunks] * x.get(i + j + half)) % FIELD_MOD;
                x.set(i + j, (even + w_odd) % FIELD_MOD);
                x.set(i + j + half, mod_sub(even, w_odd));
            }
        }
    }
}

void ntt_iter(auto& x, const arr_u128& rou_pows) {
    ntt_iter_(x, rou_pows);
}
void intt_iter(auto& x, const arr_u128& rou_pows, u128 inv_x_len) {
    ntt_iter(x, rou_pows);
    // scale by N^-1
    for (auto& el : x)
        el = (el * inv_x_len) % FIELD_MOD;
}

// Computes the recursive ntt
void ntt_recursive(span_u128 x, u128 nth_root_unity) {
    bit_reverse(x);
    ntt_recursive_(x, nth_root_unity);
}

void intt_recursive(span_u128 x, u128 nth_root_unity, u128 inv_poly_size) {
    ntt_recursive(x, nth_root_unity);
    // scale by N^-1
    for (auto& el : x)
        el = (el * inv_poly_size) % FIELD_MOD;
}