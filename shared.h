#pragma once

#include <iostream>
#include <vector>
#include <string>
#include <algorithm>
#include <chrono>
#include <gmpxx.h>

#include "serialisation.hpp"
#include "/Users/jw/Projects/mcl/include/mcl/bn.hpp"

using namespace mcl::bn;

#ifndef PAIRING_OFF
#   define PAIRING_ON
#endif

#ifdef TIMING_ON
#  define TIMING(x) x
#else
#  define TIMING(x)
#endif

typedef mpz_class bigz;
typedef long int u32; // FIXME - why signed? Rationalise usage across code
// u32 N_DECOMP = 2;
u32 N_DECOMP = 7;

// for N=4096
constexpr size_t N_ = 4096;
bigz NTH_ROU("4158865282786404163413953114870269622875596290766033564087307867933865333818");
bigz TWO_ROU("197302210312744933010843010704445784068657690384188106020011018676818793232");

// for N=8192
// constexpr size_t N_ = 8192;
// bigz NTH_ROU("197302210312744933010843010704445784068657690384188106020011018676818793232");
// bigz TWO_ROU("20619701001583904760601357484951574588621083236087856586626117568842480512645");

bigz FIELD_MODULUS("\
21888242871839275222246405745257275088548364400416034343698204186575808495617\
");

typedef struct times_and_counts {
    int calls_ntt = 0;
    int calls_intt = 0;
    int calls_ntt1 = 0;
    int calls_intt1 = 0;
    int calls_conv_to_nega = 0;
    int calls_msm = 0;
    int iter_ = 0;
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
    std::chrono::duration<double, std::milli> msm{};
    std::chrono::duration<double, std::milli> msm1{};
} times_and_counts;

// NOTE inline keyword for structs allows the struct to be used in multiple
// translation units without causing linker errors
inline times_and_counts timing = { 0 };

G1 gen1;
G2 gen2;
GT genT;


size_t N_POLYS_IN_RLWE = 2;
using matrix_double = std::vector<std::vector<mpf_class>>;
using vector_double = std::vector<mpf_class>;
using vector_bigz = std::vector<bigz>;

bigz mod_(bigz& val, const bigz& q) {
    if (val < 0 || val >= FIELD_MODULUS)
        val %= q;
    if (val < 0) {
        val += q;
    }
    return val;
}
// TODO check correctness
bigz mod_(bigz&& val, const bigz& q) {
    return mod_(val, q);
    // val = mod_(val, q);
    // return val;
}

std::string i128str(__uint128_t n) {
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
vector_bigz scalar_vec_mult(const bigz& scalar, const vector_bigz& vec, const bigz& q) {
    vector_bigz result(vec.size());
    for (size_t i = 0; i < vec.size(); i++) {
        // copy vec
        result[i] = scalar * bigz(vec[i]);
        result[i] = mod_(result[i], q);
    }
    return result;
};

void print_vector_i128(const std::vector<__uint128_t>& vec) {
    for (const auto& val : vec) {
        std::cout << i128str(val) << ", ";
    }
    std::cout << "\n";
}

void print_vector_double_old(const std::vector<double>& vec) {
    for (const auto& val : vec) {
        std::cout << val << ", ";
    }
    std::cout << "\n";
}

std::string mpf_str(mpf_class m) {
    int size = 10000;
    char* buf = new char[size];
    buf[size - 1] = '\0';
    int ret = gmp_sprintf(buf, "%Ff", m.get_mpf_t());
    if (ret > size - 1 || ret < 0) {
        throw std::runtime_error("Buffer overflow in mpf_str");
    }
    std::string s(buf);
    delete[] buf;
    return s;
}

void print_vector_double(const vector_double& vec) {
    for (const auto& val : vec) {
        std::cout << mpf_str(val) << ", ";
    }
    std::cout << "\n";
}

std::string print_to_string_mpz(bigz m) {
    int size = 10000;
    char* buf = new char[size];
    buf[size - 1] = '\0';
    int ret = gmp_sprintf(buf, "%Zd", m.get_mpz_t());
    if (ret > size - 1 || ret < 0) {
        throw std::runtime_error("Buffer overflow in print_to_string_mpz");
    }
    std::string s(buf);
    delete[] buf;
    return s;
}

void print_vector_mpz(const vector_bigz& vec) {
    for (const auto& val : vec) {
        std::cout << print_to_string_mpz(val) << ", ";
    }
    std::cout << "\n";
}
// TODO types and move somewhere
// Returns floor(log2(x))
double log2_mpz(const bigz& x) {
    if (x == 0) return -1; // or throw/handle as needed
    return static_cast<double>(mpz_sizeinbase(x.get_mpz_t(), 2) - 1);
}

mpf_class mpf_round(const mpf_class &x) {
    const mpf_class half("-0.5");
    if (x >= -1)
        return floor(x + half);   // correct for non-negative values
    else
        return ceil(x - half);    // correct for negative values
}

G1 pow_(const G1& base, const bigz& power) {
    assert(power >= 0);
    // assert(power < FIELD_MODULUS);
    // static Fr power1;
    Fr power1;
    power1.clear();
    if (power > FIELD_MODULUS) {
        std::cout << "pow_: power > FIELD_MODULUS\n";
        mpz_to_Fr(power1, power % FIELD_MODULUS);
    }
    else
        mpz_to_Fr(power1, power);
    return base * power1;
}

G2 pow2_(const G2& base, const bigz& power) {
    assert(power >= 0);
    // assert(power < FIELD_MODULUS);
    // static Fr power1;
    Fr power1;
    power1.clear();
    if (power > FIELD_MODULUS) {
        std::cout << "pow_2: power > FIELD_MODULUS\n";
        mpz_to_Fr(power1, power % FIELD_MODULUS);
    }
    else
        mpz_to_Fr(power1, power);
    return base * power1;
}

GT pow_t(const GT& base, const bigz& power) {
    assert(power >= 0);
    // assert(power < FIELD_MODULUS);
    // static Fr power1;
    Fr power1;
    power1.clear();
    if (power > FIELD_MODULUS) {
        std::cout << "pow_t: power > FIELD_MODULUS\n";
        mpz_to_Fr(power1, power % FIELD_MODULUS);
    }
    else
        mpz_to_Fr(power1, power);
    GT ret;
    GT::pow(ret, base, power1);
    return ret;
}
