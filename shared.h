#pragma once

#include <vector>
#include <string>
#include <algorithm>

// FIXME check if this is being used correctly, i.e. are some negs required?
// typedef __uint128_t i128;
typedef __int128_t i128;
constexpr i128 GROUP_MODULUS = 17180000273;
constexpr i128 FIELD_MODULUS = 1073750017;
constexpr i128 GENERATOR = 65536;

// constexpr i128 GROUP_MODULUS = 23;
// constexpr i128 FIELD_MODULUS = 11;
// constexpr i128 GENERATOR = 4;
// constexpr i128 GROUP_MODULUS = 103;
// constexpr i128 FIELD_MODULUS = 17;
// constexpr i128 GENERATOR = 64;
// constexpr i128 GROUP_MODULUS = 540431955285196831;
// constexpr unsigned int FIELD_MODULUS = 18014398509506561;
// constexpr i128 FIELD_MODULUS = 18014398509506561;
// constexpr i128 GENERATOR = 1073741824;
constexpr i128 N_POLYS_IN_RLWE = 2;
using matrix_double = std::vector<std::vector<double>>;
using vector_double = std::vector<double>;
using vector_i128 = std::vector<i128>;

// FIXME make types __uint128 so that regular modding works
i128 mod_(i128 val, i128 q) {
    val %= q;
    if (val < 0) {
        val += q;
    }
    return val;
}

vector_i128 mod_(const vector_i128& vals, i128 q) {
    vector_i128 res(vals.size());
    for (size_t i = 0; i < vals.size(); ++i) {
        res[i] = mod_(vals[i], q);
    }
    return res;
}

std::string print_to_string_i128(i128 n) {
    if (n == 0) {
        return "0";
    }
    bool neg = false;
    // if (n < 0) {
    //     neg = true;
    //     n = -n;
    // }
    std::string buf;
    while (n > 0) {
        buf += '0' + (n % 10);
        n /= 10;
    }
    if (neg) buf += '-';
    std::reverse(buf.begin(), buf.end());
    return buf;
}
vector_i128 scalar_vec_mult (i128 scalar, const vector_i128& vec, i128 q) {
    vector_i128 result(vec.size());
    for (size_t i = 0; i < vec.size(); ++i) {
        result[i] = mod_(scalar * vec[i], q);
    }
    return result;
};
