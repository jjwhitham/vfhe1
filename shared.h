#pragma once

#include <iostream>
#include <vector>
#include <string>
#include <algorithm>
#include <chrono>
// #include "gmpxx.h"
typedef __uint128_t i128;
typedef __uint32_t u32;
constexpr u32 N_DECOMP = 2;

#ifdef TIMING_ON
#  define TIMING(x) x
#else
#  define TIMING(x)
#endif

typedef struct {
    int calls_ntt = 0;
    int calls_intt = 0;
    int calls_ntt1 = 0;
    int calls_intt1 = 0;
    int calls_conv_to_nega = 0;
    int calls_get_hash_sec = 0;
    i128 iter_ = 0;
    std::chrono::duration<double, std::milli> verify{};
    std::chrono::duration<double, std::milli> proof{};
    std::chrono::duration<double, std::milli> controller{};
    std::chrono::duration<double, std::milli> plant{};
    std::chrono::duration<double, std::milli> total{};
    std::chrono::duration<double, std::milli> loop{};
    std::chrono::duration<double, std::milli> ntt{};
    std::chrono::duration<double, std::milli> intt{};
    std::chrono::duration<double, std::milli> ntt1{};
    std::chrono::duration<double, std::milli> intt1{};
    std::chrono::duration<double, std::milli> get_hash_sec{};
} times_and_counts;

// NOTE inline keyword for structs allows the struct to be used in multiple
// translation units without causing linker errors
inline times_and_counts timing = { 0 };

// constexpr i128 N_ = 2;
// constexpr i128 GROUP_MODULUS = 11; // p
// constexpr i128 FIELD_MODULUS = 5; // q
// constexpr i128 GENERATOR = 4; // g

// q_pow = 30.000011008191272;
// constexpr i128 N_ = 4096;
// constexpr i128 GROUP_MODULUS = 17180000273;
// constexpr i128 FIELD_MODULUS = 1073750017;
// constexpr i128 GENERATOR = 65536;
// constexpr i128 NTH_ROU = 625534531;
// constexpr i128 TWO_ROU = 996876704;

// TODO check N/q
// N = 2^4, q > 2^4
// constexpr i128 GROUP_MODULUS = 103;
// constexpr i128 FIELD_MODULUS = 17;
// constexpr i128 GENERATOR = 64;

// N = 2^12, q > 2^54
constexpr i128 N_ = 4096;
constexpr i128 GROUP_MODULUS = 540431955285196831;
constexpr i128 FIELD_MODULUS = 18014398509506561;
constexpr i128 GENERATOR = 1073741824;
constexpr i128 NTH_ROU = 5194839201355896;
constexpr i128 TWO_ROU = 9455140237568613;

/* constexpr mpz GROUP_MODULUS = \
4893786199427765557762591497536453196573320618411278635472176938119840441271\
2502191283384375987166735470567071036504573974100823080776092800473413359198247\
6239041963014833898372354040701042388893258698981186938413740280916319621417623\
7471114636814326740929854005470025012660047589916149139703043020351672528846339\
84145576703811773398576717571852999
constexpr mpz FIELD_MODULUS = 340282366920938463463374607431767867393
constexpr mpz GENERATOR = \
4817716213549514035336014710666703183713025673847494424413528199993357208789\
8913394575489847649253866038828334491235671772417871456974026431473237742550247\
9497578100053448330287050686965459146273651762955474911950522760806904335945868\
4003625854903563191461961625964949673337592873574451932801584164249941551882456\
08741421710111223106441384922396297 */

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

std::string i128str(i128 n) {
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
        std::cout << i128str(val) << ", ";
    }
    std::cout << "\n";
}

void print_vector_double(const vector_double& vec) {
    for (const auto& val : vec) {
        std::cout << val << ", ";
    }
    std::cout << "\n";
}