#include <vector>
#include <iostream>
#include "atcoder/convolution.hpp" // Download from AtCoder ACL

using namespace std;
using namespace atcoder;
typedef __int128_t i128;

// const i128 MOD = 998244353;
// NOTE Must be a suitable NTT modulus (m * 2^k + 1)
const i128 MOD = 1 * (1 << 4) + 1; // 17

std::string print_to_string_i128(i128 n) {
    if (n == 0) {
        return "0";
    }
    bool neg = false;
    if (n < 0) {
        neg = true;
        n = -n;
    }
    std::string buf;
    while (n > 0) {
        buf += '0' + (n % 10);
        n /= 10;
    }
    if (neg) buf += '-';
    std::reverse(buf.begin(), buf.end());
    return buf;
}

// Negacyclic convolution: c(x) = a(x) * b(x) mod (x^n + 1)
vector<i128> negacyclic_convolution(const vector<i128>& a, const vector<i128>& b) {
    i128 n = a.size();
    vector<i128> a_pad = a, b_pad = b;
    a_pad.resize(2 * n);
    b_pad.resize(2 * n);
    vector<i128> conv = convolution<MOD>(a_pad, b_pad);
    for (i128 c : conv) std::cout << print_to_string_i128(c) << " ";
    std::cout << "\n";
    vector<i128> res(n);
    for (i128 i = 0; i < n; ++i) {
        // (conv[i] - conv[i + n]) mod MOD
        res[i] = (conv[i] - conv[i + n] + MOD) % MOD;
    }
    return res;
}

// class arr {
// private:
//     int size_;
//     int* arr_;
// public:
//     arr() : size_(0) {}
//     arr(size_t size) : size_(size) {
//         arr_ = new int[size]();
//     }
//     ~arr() = default;
// };

int main() {
    vector<i128> a = {1, 2};
    vector<i128> b = {3, 4};
    // vector<i128> b = {5, 6, 7, 8};
    vector<i128> c = convolution<MOD>(a, b);
    vector<i128> c1 = negacyclic_convolution(a, b);
    cout << "Convolution result: ";
    for (i128 x : c) cout << print_to_string_i128(x) << " ";
    cout << "\n";
    cout << "negacyclic convolution result: ";
    for (i128 x : c1) cout << print_to_string_i128(x) << " ";
    cout << "\n";
    // arr ar(2);
    // arr aconv = convolution<MOD>(a, b);
}

// TODO change vector to array, or modify enc_types to take vector?
// TODO add isNTT flag
// TODO add hashed
// TODO add convolution
// TODO