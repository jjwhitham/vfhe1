#pragma once

#include <iostream>
#include <vector>
#include <string>
#include <algorithm>
// #include "gmpxx.h"

typedef __uint128_t i128;
typedef __uint32_t u32;
constexpr u32 N_DECOMP = 1;

// constexpr i128 N_ = 2;
// constexpr i128 GROUP_MODULUS = 11; // p
// constexpr i128 FIELD_MODULUS = 5; // q
// constexpr i128 GENERATOR = 4; // g



// #include <gmpxx.h>

// typedef __uint128_t u128;

// void find_q(mpz_class& q, mpz_class const& N, mpz_class const& q_init) {
//     int i = 0;
//     mpz_class k = (q_init - 1) / (2 * N) + 1;
//     q = q_init;
//     if (mpz_probab_prime_p(q.get_mpz_t(), 25) && (q % (2 * N) == 1))
//         return;
//     while (i < 100) {
//         mpz_nextprime(q.get_mpz_t(), q.get_mpz_t());
//         if (q % (2 * N) == 1)
//             return;
//         k += 1;
//         // q = 1 + 2 * N * k; // FIXME drop the 2?
//         i += 1;
//     }
// }

// q_pow = 30.000000044339277;
// constexpr i128 N_ = 16;
// constexpr i128 GROUP_MODULUS = 36507223139;
// constexpr i128 FIELD_MODULUS = 1073741857;
// constexpr i128 GENERATOR = 17179869184;

// // q_pow = 30.000011008191272;
// constexpr i128 N_ = 256;
// constexpr i128 GROUP_MODULUS = 17180000273;
// constexpr i128 FIELD_MODULUS = 1073750017;
// constexpr i128 GENERATOR = 65536;

// // q_pow = 30.000011008191272;
// constexpr i128 N_ = 1024;
// constexpr i128 GROUP_MODULUS = 17180000273;
// constexpr i128 FIELD_MODULUS = 1073750017;
// constexpr i128 GENERATOR = 65536;


// q_pow = 30.000011008191272;
// constexpr i128 N_ = 2048;
// constexpr i128 GROUP_MODULUS = 17180000273;
// constexpr i128 FIELD_MODULUS = 1073750017;
// constexpr i128 GENERATOR = 65536;

// q_pow = 30.000011008191272;
constexpr i128 N_ = 4096;
constexpr i128 GROUP_MODULUS = 17180000273;
constexpr i128 FIELD_MODULUS = 1073750017;
constexpr i128 GENERATOR = 65536;
constexpr i128 NTH_ROU = 625534531;
constexpr i128 TWO_ROU = 996876704;




// q_pow = 5;
// constexpr i128 N_ = 4096;
// constexpr i128 GROUP_MODULUS = 327689;
// constexpr i128 FIELD_MODULUS = 40961;
// constexpr i128 GENERATOR = 256;
// constexpr i128 NTH_ROU = 18088;
// constexpr i128 TWO_ROU = 243;

// // q_pow = 10;
// constexpr i128 N_ = 4096;
// constexpr i128 GROUP_MODULUS = 327689;
// constexpr i128 FIELD_MODULUS = 40961;
// constexpr i128 GENERATOR = 256;

// // q_pow = 20;
// constexpr i128 N_ = 4096;
// constexpr i128 GROUP_MODULUS = 2146307;
// constexpr i128 FIELD_MODULUS = 1073153;
// constexpr i128 GENERATOR = 4;

// q_pow = 25;
// constexpr i128 N_ = 4096;
// constexpr i128 GROUP_MODULUS = 674201621;
// constexpr i128 FIELD_MODULUS = 33710081;
// constexpr i128 GENERATOR = 1048576;

// // q_pow = 30;

// constexpr i128 N_ = 4096;
// constexpr i128 GROUP_MODULUS = 17180000273;
// constexpr i128 FIELD_MODULUS = 1073750017;
// constexpr i128 GENERATOR = 65536;


// N = 2^12, q > 2^30
// constexpr i128 GROUP_MODULUS = 17180000273;
// constexpr i128 FIELD_MODULUS = 1073750017;
// constexpr i128 GENERATOR = 65536;

// N = 2^12, q > 2^12
// constexpr i128 GROUP_MODULUS = 327689;
// constexpr i128 FIELD_MODULUS = 40961;
// constexpr i128 GENERATOR = 256;

// N = 2^3, q > 2^3
// constexpr i128 GROUP_MODULUS = 23;
// constexpr i128 FIELD_MODULUS = 11;
// constexpr i128 GENERATOR = 4;

// TODO check N/q
// N = 2^4, q > 2^4
// constexpr i128 GROUP_MODULUS = 103;
// constexpr i128 FIELD_MODULUS = 17;
// constexpr i128 GENERATOR = 64;

// N = 2^12, q > 2^54
// constexpr i128 N_ = 4096;
// constexpr i128 GROUP_MODULUS = 540431955285196831;
// constexpr i128 FIELD_MODULUS = 18014398509506561;
// constexpr i128 GENERATOR = 1073741824;
// constexpr i128 NTH_ROU = 5194839201355896;
// constexpr i128 TWO_ROU = 9455140237568613;


constexpr i128 N_POLYS_IN_RLWE = 2;
using matrix_double = std::vector<std::vector<double>>;
using vector_double = std::vector<double>;
using vector_i128 = std::vector<i128>;

// FIXME make types __uint128 so that regular modding works
i128 mod_(i128 val, i128 q) {
    val %= q;
    // if (val < 0) {
    //     val += q;
    // }
    return val;
}

// mpz_class mod_(mpz_class val, mpz_class q) {
//     val %= q;
//     if (val < 0) {
//         val += q;
//     }
//     return val;
// }


vector_i128 mod_(const vector_i128& vals, i128 q) {
    vector_i128 res(vals.size());
    for (size_t i = 0; i < vals.size(); i++) {
        res[i] = mod_(vals[i], q);
    }
    return res;
}

std::string print_to_string_i128(i128 n) {
    if (n == 0) {
        return "0";
    }
    // bool neg = false;
    // if (n < 0) {
    //     neg = true;
    //     n = -n;
    // }
    std::string buf;
    while (n > 0) {
        buf += '0' + (n % 10);
        n /= 10;
    }
    // if (neg) buf += '-';
    std::reverse(buf.begin(), buf.end());
    return buf;
}
vector_i128 scalar_vec_mult(i128 scalar, const vector_i128& vec, i128 q) {
    vector_i128 result(vec.size());
    for (size_t i = 0; i < vec.size(); i++) {
        result[i] = mod_(scalar * vec[i], q);
    }
    return result;
};

void print_vector_i128(const vector_i128& vec) {
    for (const auto& val : vec) {
        std::cout << print_to_string_i128(val) << ", ";
    }
    std::cout << "\n";
}

void print_vector_double(const vector_double& vec) {
    for (const auto& val : vec) {
        std::cout << val << ", ";
    }
    std::cout << "\n";
}