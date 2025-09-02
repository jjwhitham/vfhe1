#include <sodium.h>
#include <iostream>
#include <iomanip>
#include <cstring>

// [const char | char const][* | * const][* | * const][* | * const]
// uchar_cppc
using uchar = unsigned char;
using uchar_p = unsigned char*;
using uchar_c = const unsigned char;
using uch_cp = const unsigned char*;
using u128 = __uint128_t;
constexpr size_t N_BYTES = crypto_core_ristretto255_BYTES;

uch_cp conv_to_256(auto in) {
    uchar_p out = new uchar[N_BYTES];
    std::memset(out, 0, N_BYTES);
    std::memcpy(out, static_cast<void*>(&in), sizeof(in));
    return out;
}

static void print_hex(const unsigned char *buf, size_t len) {
    for (size_t i = 0; i < len; ++i)
        std::cout << std::hex << std::setw(2) << std::setfill('0') << (int)buf[i];
    std::cout << std::dec << "\n";
}

int main() {
    if (sodium_init() < 0) return 1;

    unsigned char a[32], b[32];
    unsigned char P[32], L[32], R[32], C[32];

    // generate random scalars (Ristretto scalar API)
    crypto_core_ristretto255_scalar_random(a);
    crypto_core_ristretto255_scalar_random(b);

    // P = b * B   (Ristretto basepoint)
    if (crypto_scalarmult_ristretto255_base(P, b) != 0) return 1;

    // L = a * P
    if (crypto_scalarmult_ristretto255(L, a, P) != 0) return 1;

    // C = a * b (scalar mul mod group order)
    crypto_core_ristretto255_scalar_mul(C, a, b);
    // R = C * B
    if (crypto_scalarmult_ristretto255_base(R, C) != 0) return 1;

    std::cout << "L: "; print_hex(L, 32);
    std::cout << "R: "; print_hex(R, 32);
    std::cout << (std::memcmp(L, R, 32) == 0 ? "OK\n" : "FAIL\n");

    std::memset(a, 0, 32);
    int x = 1234566;
    std::memcpy(a, static_cast<void*>(&x), 4);
    std::memset(b, 0, 32);
    x = 90812547;
    std::memcpy(b, static_cast<void*>(&x), 4);

    unsigned char lhs[32], rhs[32], aB[32], bB[32];
    // (a + b) * B
    crypto_core_ristretto255_scalar_add(lhs, a, b);
    crypto_scalarmult_ristretto255_base(lhs, lhs);
    // aB + bB
    crypto_scalarmult_ristretto255_base(aB, a);
    crypto_scalarmult_ristretto255_base(bB, b);
    crypto_core_ristretto255_add(rhs, aB, bB);

    std::cout << (std::memcmp(lhs, rhs, 32) == 0 ? "OK\n" : "FAIL\n");


    int xx = 42;
    uint64_t y = 12948591825;
    __uint64_t z = 1985749598178;
    uch_cp xxx = conv_to_256(xx);
    uch_cp yy = conv_to_256(y);
    uch_cp zz = conv_to_256(z);
    return 0;
}

// #include <sodium.h>
// #include <iostream>
// #include <iomanip>
// #include <cstring>
// #include <vector>
// #include <array>
// #include <chrono>

// static void print_hex(const unsigned char *buf, size_t len) {
//     for (size_t i = 0; i < len; ++i)
//         std::cout << std::hex << std::setw(2) << std::setfill('0') << (int)buf[i];
//     std::cout << std::dec << "\n";
// }

// int main() {
//     if (sodium_init() < 0) return 1;

//     constexpr size_t BUFSZ = 4096;
//     constexpr size_t WARMUP = 1024;
//     constexpr size_t ITERS  = 20000;

//     // pre-generate random scalars and points
//     std::vector<std::array<unsigned char,32>> scalars(BUFSZ);
//     std::vector<std::array<unsigned char,32>> points(BUFSZ);
//     for (size_t i = 0; i < BUFSZ; ++i) {
//         crypto_core_ristretto255_scalar_random(scalars[i].data());
//         crypto_scalarmult_ristretto255_base(points[i].data(), scalars[i].data());
//     }

//     volatile unsigned char checksum = 0;

//     // warmup arbitrary scalar * point
//     {
//         unsigned char out[32];
//         for (size_t i = 0; i < WARMUP; ++i) {
//             size_t idx = i % BUFSZ;
//             crypto_scalarmult_ristretto255(out, scalars[idx].data(), points[(idx+1)%BUFSZ].data());
//             checksum ^= out[0];
//         }
//     }

//     // measure scalar * point
//     double secs1;
//     {
//         unsigned char out[32];
//         auto t0 = std::chrono::steady_clock::now();
//         for (size_t i = 0; i < ITERS; ++i) {
//             size_t idx = i % BUFSZ;
//             crypto_scalarmult_ristretto255(out, scalars[idx].data(), points[(idx+1)%BUFSZ].data());
//             checksum ^= out[0];
//         }
//         auto t1 = std::chrono::steady_clock::now();
//         secs1 = std::chrono::duration<double>(t1 - t0).count();
//         double ops_per_sec = double(ITERS) / secs1;
//         std::cout << "scalar * point (crypto_scalarmult_ristretto255):\n";
//         std::cout << "  iterations = " << ITERS << ", elapsed = " << secs1 << " s\n";
//         std::cout << "  ops/sec = " << std::fixed << std::setprecision(0) << ops_per_sec
//                   << "  ( " << std::setprecision(3) << (1e6 / ops_per_sec) << " us/op )\n\n";
//     }

//     // warmup scalar * base
//     {
//         unsigned char out[32];
//         for (size_t i = 0; i < WARMUP; ++i) {
//             size_t idx = i % BUFSZ;
//             crypto_scalarmult_ristretto255_base(out, scalars[idx].data());
//             checksum ^= out[0];
//         }
//     }

//     // measure scalar * base
//     double secs2;
//     {
//         unsigned char out[32];
//         auto t0 = std::chrono::steady_clock::now();
//         for (size_t i = 0; i < ITERS; ++i) {
//             size_t idx = i % BUFSZ;
//             crypto_scalarmult_ristretto255_base(out, scalars[idx].data());
//             checksum ^= out[0];
//         }
//         auto t1 = std::chrono::steady_clock::now();
//         secs2 = std::chrono::duration<double>(t1 - t0).count();
//         double ops_per_sec = double(ITERS) / secs2;
//         std::cout << "scalar * base (crypto_scalarmult_ristretto255_base):\n";
//         std::cout << "  iterations = " << ITERS << ", elapsed = " << secs2 << " s\n";
//         std::cout << "  ops/sec = " << std::fixed << std::setprecision(0) << ops_per_sec
//                   << "  ( " << std::setprecision(3) << (1e6 / ops_per_sec) << " us/op )\n\n";
//     }

//     // warmup chained exponentiation (a*(b*B))
//     {
//         unsigned char p[32], r[32];
//         for (size_t i = 0; i < WARMUP; ++i) {
//             size_t idx = i % BUFSZ;
//             crypto_scalarmult_ristretto255_base(p, scalars[idx].data());
//             crypto_scalarmult_ristretto255(r, scalars[(idx+1)%BUFSZ].data(), p);
//             checksum ^= r[0];
//         }
//     }

//     // measure chained exponentiation: one "exp" = two scalar multiplications (base then mul)
//     double secs3;
//     {
//         unsigned char p[32], r[32];
//         auto t0 = std::chrono::steady_clock::now();
//         for (size_t i = 0; i < ITERS; ++i) {
//             size_t idx = i % BUFSZ;
//             crypto_scalarmult_ristretto255_base(p, scalars[idx].data());         // P = b*B
//             crypto_scalarmult_ristretto255(r, scalars[(idx+1)%BUFSZ].data(), p); // R = a*P
//             checksum ^= r[0];
//         }
//         auto t1 = std::chrono::steady_clock::now();
//         secs3 = std::chrono::duration<double>(t1 - t0).count();
//         double exps_per_sec = double(ITERS) / secs3;
//         std::cout << "chained exponentiation (a*(b*B)) [2 scalar-mults per exp]:\n";
//         std::cout << "  iterations = " << ITERS << ", elapsed = " << secs3 << " s\n";
//         std::cout << "  exp/sec = " << std::fixed << std::setprecision(0) << exps_per_sec
//                   << "  ( " << std::setprecision(3) << (1e6 / exps_per_sec) << " us/exp )\n";
//         std::cout << "  scalar-mults/sec (raw) = " << std::fixed << std::setprecision(0) << (exps_per_sec * 2) << "\n\n";
//     }

//     // print checksum to avoid optimization removal
//     std::cout << "checksum: " << std::hex << int(checksum) << std::dec << "\n";
//     return 0;
// }