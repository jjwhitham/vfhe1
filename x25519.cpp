#include <sodium.h>
#include <iostream>
#include <iomanip>
#include <cstring>

static void print_hex(const unsigned char *buf, size_t len) {
    for (size_t i = 0; i < len; ++i) std::cout << std::hex << std::setw(2) << std::setfill('0') << (int)buf[i];
    std::cout << std::dec << "\n";
}

int main() {
    if (sodium_init() < 0) return 1;

    // X25519 keys are 32 bytes (256 bits)
    unsigned char alice_pk[32], alice_sk[32];
    unsigned char bob_pk[32], bob_sk[32];

    // Generate keypairs (libsodium will produce proper X25519 keys)
    crypto_kx_keypair(alice_pk, alice_sk);
    crypto_kx_keypair(bob_pk, bob_sk);

    std::cout << "alice_sk: "; print_hex(alice_sk, 32);
    std::cout << "alice_pk: "; print_hex(alice_pk, 32);
    std::cout << "bob_sk:   "; print_hex(bob_sk, 32);
    std::cout << "bob_pk:   "; print_hex(bob_pk, 32);

    // Group operation: scalar multiplication (Diffie-Hellman)
    unsigned char shared1[32], shared2[32];
    if (crypto_scalarmult(shared1, alice_sk, bob_pk) != 0) return 1;
    if (crypto_scalarmult(shared2, bob_sk, alice_pk) != 0) return 1;

    std::cout << "shared1:  "; print_hex(shared1, 32);
    std::cout << "shared2:  "; print_hex(shared2, 32);
    std::cout << "shared match: " << (std::memcmp(shared1, shared2, 32) == 0) << "\n";

    // Exponentiation example: compute R = k2 * (k1 * base)
    unsigned char k1[32], k2[32];
    randombytes_buf(k1, 32);
    randombytes_buf(k2, 32);

    unsigned char P[32], R[32];
    // P = k1 * basepoint
    if (crypto_scalarmult_base(P, k1) != 0) return 1;
    // R = k2 * P  -> (k2 * k1) * base
    if (crypto_scalarmult(R, k2, P) != 0) return 1;

    std::cout << "k1:       "; print_hex(k1, 32);
    std::cout << "k2:       "; print_hex(k2, 32);
    std::cout << "P = k1*B: "; print_hex(P, 32);
    std::cout << "R = k2*P: "; print_hex(R, 32);

    // Verify associativity via scalar multiplication with product scalar:
    // compute kprod = k1 * k2 (mod group order) is not directly meaningful here;
    // but check kprod*base == R by repeated scalar-mult: compute S = (kprod)*base using libsodium isn't provided,
    // so we demonstrate composition via scalarmult chaining above.

    return 0;
}