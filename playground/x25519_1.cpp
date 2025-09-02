// ...existing code...
#include <sodium.h>
#include <iostream>
#include <iomanip>
#include <cstring>
#include <vector>
#include <chrono>

static void print_hex(const unsigned char *buf, size_t len) {
    for (size_t i = 0; i < len; ++i) std::cout << std::hex << std::setw(2) << std::setfill('0') << (int)buf[i];
    std::cout << std::dec << "\n";
}

int main() {
    if (sodium_init() < 0) return 1;

    // Warmup / benchmark params
    constexpr size_t BUFSZ = 4096;      // number of pre-generated inputs to cycle through
    constexpr size_t WARMUP = 1024;
    constexpr size_t ITERS  = 20000;    // total measured iterations

    // Pre-generate random scalars and points to avoid RNG overhead in measured loop
    std::vector<std::array<unsigned char,32>> scalars(BUFSZ);
    std::vector<std::array<unsigned char,32>> points(BUFSZ);
    for (size_t i = 0; i < BUFSZ; ++i) {
        randombytes_buf(scalars[i].data(), 32);
        // create a valid point by multiplying base by a random scalar
        crypto_scalarmult_base(points[i].data(), scalars[i].data());
    }

    // checksum to prevent optimizer from removing calls
    volatile unsigned char check = 0;

    // Warmup: scalar mult (scalar * arbitrary point)
    {
        unsigned char out[32];
        for (size_t i = 0; i < WARMUP; ++i) {
            size_t idx = i % BUFSZ;
            crypto_scalarmult(out, scalars[idx].data(), points[(idx+1)%BUFSZ].data());
            check ^= out[0];
        }
    }

    // Measured: crypto_scalarmult (scalar * arbitrary point)
    double secs = 0;
    {
        unsigned char out[32];
        auto t0 = std::chrono::steady_clock::now();
        for (size_t i = 0; i < ITERS; ++i) {
            size_t idx = i % BUFSZ;
            crypto_scalarmult(out, scalars[idx].data(), points[(idx+1)%BUFSZ].data());
            check ^= out[0];
        }
        auto t1 = std::chrono::steady_clock::now();
        secs = std::chrono::duration<double>(t1 - t0).count();
        double ops_per_sec = double(ITERS) / secs;
        std::cout << "crypto_scalarmult (scalar * point):\n";
        std::cout << "  iterations = " << ITERS << ", elapsed = " << secs << " s\n";
        std::cout << "  ops/sec = " << std::fixed << std::setprecision(0) << ops_per_sec
                  << "  ( " << std::setprecision(3) << (1e6 / ops_per_sec) << " us/op )\n\n";
    }

    // Warmup: scalar mult base
    {
        unsigned char out[32];
        for (size_t i = 0; i < WARMUP; ++i) {
            size_t idx = i % BUFSZ;
            crypto_scalarmult_base(out, scalars[idx].data());
            check ^= out[0];
        }
    }

    // Measured: crypto_scalarmult_base (scalar * base)
    {
        unsigned char out[32];
        auto t0 = std::chrono::steady_clock::now();
        for (size_t i = 0; i < ITERS; ++i) {
            size_t idx = i % BUFSZ;
            crypto_scalarmult_base(out, scalars[idx].data());
            check ^= out[0];
        }
        auto t1 = std::chrono::steady_clock::now();
        double secs2 = std::chrono::duration<double>(t1 - t0).count();
        double ops_per_sec = double(ITERS) / secs2;
        std::cout << "crypto_scalarmult_base (scalar * base):\n";
        std::cout << "  iterations = " << ITERS << ", elapsed = " << secs2 << " s\n";
        std::cout << "  ops/sec = " << std::fixed << std::setprecision(0) << ops_per_sec
                  << "  ( " << std::setprecision(3) << (1e6 / ops_per_sec) << " us/op )\n\n";
    }

    // Demonstrate chained multiplication (k2 * (k1 * base)) throughput:
    // warmup
    {
        unsigned char p[32], r[32];
        for (size_t i = 0; i < WARMUP; ++i) {
            size_t idx = i % BUFSZ;
            crypto_scalarmult_base(p, scalars[idx].data());       // P = k1*B
            crypto_scalarmult(r, scalars[(idx+1)%BUFSZ].data(), p); // R = k2*P
            check ^= r[0];
        }
    }

    // measure chained exponentiation (two scalar multiplies per "exp")
    {
        unsigned char p[32], r[32];
        auto t0 = std::chrono::steady_clock::now();
        for (size_t i = 0; i < ITERS; ++i) {
            size_t idx = i % BUFSZ;
            crypto_scalarmult_base(p, scalars[idx].data());         // P = k1*B
            crypto_scalarmult(r, scalars[(idx+1)%BUFSZ].data(), p); // R = k2*P
            check ^= r[0];
        }
        auto t1 = std::chrono::steady_clock::now();
        double secs3 = std::chrono::duration<double>(t1 - t0).count();
        double ops_per_sec = double(ITERS) / secs3; // counts chained ops as 1 "exp"
        std::cout << "chained exponentiation (k2*(k1*B)) [two scalar mults per exp]:\n";
        std::cout << "  iterations = " << ITERS << ", elapsed = " << secs3 << " s\n";
        std::cout << "  exp/sec = " << std::fixed << std::setprecision(0) << ops_per_sec
                  << "  ( " << std::setprecision(3) << (1e6 / ops_per_sec) << " us/exp )\n";
        // also report raw scalar-mults/sec (2 * exp/sec)
        std::cout << "  scalar-mults/sec (raw) = " << std::fixed << std::setprecision(0) << (ops_per_sec * 2) << "\n\n";
    }

    // print checksum to avoid optimizer removing calls
    std::cout << "checksum: " << std::hex << int(check) << std::dec << "\n";

    return 0;
}
// ...existing code...