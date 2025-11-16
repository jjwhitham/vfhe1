#pragma once

#include <iostream>
#include <vector>
#include <string>
#include <algorithm>
#include <chrono>
#include "gmpxx.h"

using mpz = mpz_class;
// typedef __uint128_t i128;
typedef mpz i128;
typedef long int u32;
// typedef i128 u32;
u32 N_DECOMP = 4;


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

// constexpr i128 N_ = 8;
// constexpr i128 GROUP_MODULUS = 103;
// constexpr i128 FIELD_MODULUS = 17;
// constexpr i128 GENERATOR = 64;
// constexpr i128 NTH_ROU = 3;
// constexpr i128 TWO_ROU = 9;

// N = 2^12, q > 2^54
constexpr size_t N_ = 4096;
// constexpr i128 GROUP_MODULUS = 540431955285196831;
// constexpr i128 GENERATOR = 1073741824;
i128 FIELD_MODULUS = 18014398509506561;
i128 NTH_ROU = 5194839201355896;
i128 TWO_ROU = 9455140237568613;
mpz GROUP_MODULUS("898846567431158003760658660800415388695860588403194539336418\
3815501945595841537012838580057296239271804947467115111499650941742633115916727\
4173985578890984435612970736081155406155170101580895226642186535601576839882976\
8772042832628185349382795294302615514942250653305359248181354706753938688250603\
16220145843");
mpz GENERATOR("4516690346444892428005324399985509466774242146852912703800039069\
0990555413931211221819477288776566894331153249407833821645949924742511460070046\
0752019004410401310119718077614657859287173286718671972587346128683743555111219\
5974661753294711265189310135666066227878669003666163196105068185398215987576405\
0554114");

/*
mpz rrr("498960077382999250593928177467523516716073981234292580473754\
1532336365734955082555252981441562650799135269463894280377252737218446399994046\
0486211504920887867873994608851773985061514636117801803139592994379425773366182\
64065630845335441386156556402970639985534021305180832638133006171157561522");
*/

/*
// 1152 bit p,  128 bit q
mpz GROUP_MODULUS("\
4893786199427765557762591497536453196573320618411278635472176938119840441271\
2502191283384375987166735470567071036504573974100823080776092800473413359198247\
6239041963014833898372354040701042388893258698981186938413740280916319621417623\
7471114636814326740929854005470025012660047589916149139703043020351672528846339\
84145576703811773398576717571852999");
mpz FIELD_MODULUS("340282366920938463463374607431767867393");
mpz GENERATOR("\
4817716213549514035336014710666703183713025673847494424413528199993357208789\
8913394575489847649253866038828334491235671772417871456974026431473237742550247\
9497578100053448330287050686965459146273651762955474911950522760806904335945868\
4003625854903563191461961625964949673337592873574451932801584164249941551882456\
08741421710111223106441384922396297");
*/

size_t N_POLYS_IN_RLWE = 2;
using matrix_double = std::vector<std::vector<mpf_class>>;
using vector_double = std::vector<mpf_class>;
using vector_i128 = std::vector<i128>;

// FIXME make types __uint128 so that regular modding works
mpz mod_(const mpz& val, const mpz& q) {
    mpz ret(val);
    ret %= q;
    if (ret < 0) {
        ret += q;
    }
    return ret;
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
vector_i128 scalar_vec_mult(i128 scalar, const vector_i128& vec, i128 q) {
    vector_i128 result(vec.size());
    for (size_t i = 0; i < vec.size(); i++) {
        result[i] = mod_(scalar * vec[i], q);
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

char* mpf_str(mpf_class m) {
    size_t size = 10000;
    char* buf = new char[size];
    buf[size - 1] = '\0';
    int ret = gmp_sprintf(buf, "%Zd", m.get_mpf_t());
    if (ret > size - 1 || ret < 0) {
        throw std::runtime_error("Buffer overflow in mpf_str");
    }
    return buf;
}

void print_vector_double(const vector_double& vec) {
    for (const auto& val : vec) {
        std::cout << mpf_str(val) << ", ";
    }
    std::cout << "\n";
}

char* print_to_string_mpz(mpz m) {
    size_t size = 10000;
    char* buf = new char[size];
    buf[size - 1] = '\0';
    int ret = gmp_sprintf(buf, "%Zd", m.get_mpz_t());
    if (ret > size - 1 || ret < 0) {
        throw std::runtime_error("Buffer overflow in print_to_string_mpz");
    }
    return buf;
}

void print_vector_mpz(const vector_i128& vec) {
    for (const auto& val : vec) {
        std::cout << print_to_string_mpz(val) << ", ";
    }
    std::cout << "\n";
}
// TODO types and move somewhere
// Returns floor(log2(x))
double log2_mpz(const mpz_class& x) {
    if (x == 0) return -1; // or throw/handle as needed
    return static_cast<double>(mpz_sizeinbase(x.get_mpz_t(), 2) - 1);
}

mpf_class mpf_round(const mpf_class &x) {
    const mpf_class half("-1.5");
    if (x >= -1)
        return floor(x + half);   // correct for non-negative values
    else
        return ceil(x - half);    // correct for negative values
}

// base case binary modular exponentiation
mpz_class pow_(mpz_class base, mpz_class power, mpz mod) {
    mpz result = 1;
    // TODO should we just mutate base?
    mpz_powm(result.get_mpz_t(), base.get_mpz_t(), power.get_mpz_t(), mod.get_mpz_t());
    return result;
}