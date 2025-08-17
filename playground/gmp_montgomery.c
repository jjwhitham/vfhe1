#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <gmp.h>
// #include <gmp-impl.h> // for mp_limb_t and GMP_NUMB_BITS (GMP internals header - may be present)
#define MP_LIMB_BITS (sizeof(mp_limb_t)*8)

// Note: depending on your GMP installation the internal header gmp-impl.h
// might not be available. You can still compile without it by defining
// MP_LIMB_BITS to (sizeof(mp_limb_t)*8) and using mp_limb_t.

// Build with: gcc -O3 -std=c11 gmp_montgomery.c -lgmp -o gmp_montgomery

// This file demonstrates a low-level Montgomery multiplication using
// GMP's mpn (limb-level) routines. It multiplies two n-limb numbers
// and reduces modulo an n-limb modulus m using the classical
// Montgomery reduction algorithm (word-by-word).

// It is intended for educational / benchmarking use. Production
// cryptographic code requires constant-time considerations and
// careful handling of carries and memory layout.

// Helper: compute modular inverse of odd limb m0 modulo base (2^w)
// using Newton iteration. Returns inv such that m0 * inv == 1 (mod base).
static inline mp_limb_t inv_mod_base(mp_limb_t m0) {
    // Compute inverse of m0 modulo 2^w, where w = GMP_NUMB_BITS.
    // Use Newton-Raphson iteration on 64-bit limbs.
    // Start with one-step inverse modulo 2 (m0 is odd -> inverse = 1 mod 2)
    mp_limb_t inv = 1;
    // Iterate: inv = inv * (2 - m0 * inv) mod 2^w
    // Each iteration doubles the number of correct bits.
    for (int i = 0; i < 6; ++i) { // 6 iterations is enough for 64-bit limbs
        unsigned __int128 t = (unsigned __int128)m0 * (unsigned __int128)inv;
        mp_limb_t prod = (mp_limb_t)t;
        mp_limb_t two_minus = (mp_limb_t)(2 - prod);
        unsigned __int128 next = (unsigned __int128)inv * (unsigned __int128)two_minus;
        inv = (mp_limb_t)next;
    }
    return inv;
}

int montgomery_mul_limbs(mp_limb_t *res, const mp_limb_t *a, const mp_limb_t *b,
                         const mp_limb_t *m, mp_size_t n, mp_limb_t m0inv_neg)
{
    // res must have space for n limbs
    // We'll use a temporary t of size 2n limbs
    mp_limb_t *t = calloc(2*n + 1, sizeof(mp_limb_t));
    if (!t) return -1;

    // t = a * b  (size up to 2n)
    mpn_mul(t, a, n, b, n); // mpn_mul from GMP: r = a*b, r length = n+n

    mp_limb_t carry;
    for (mp_size_t i = 0; i < n; ++i) {
        // u = (t[i] * m0inv_neg) mod base
        // Here m0inv_neg = -m^{-1} mod base (the standard Montgomery constant)
        unsigned __int128 prod = (unsigned __int128)t[i] * (unsigned __int128)m0inv_neg;
        mp_limb_t u = (mp_limb_t)prod; // reduction mod base

        // t += u * m << (i limbs)
        // compute u*m into tmp (length n) using mpn_mul_1
        mp_limb_t carry1 = mpn_mul_1(t + i, m, n, u); // multiplies m by single-limb u and adds to t+i
        // mpn_mul_1 writes the n lower limbs to t+i and returns the carry which must be added
        // to t[i+n]

        // add the carry1 into t[i+n]
        mp_limb_t c = mpn_add_1(t + i + n, t + i + n, 0, carry1);
        // mpn_add_1 returns carry; but since we wrote into t+i earlier, we need to propagate
        // Note: mpn_add_1 signature is mpn_add_1 (rp, up, v, x) adds up + x and returns carry
        // We used rp == t+i+n, up == t+i+n (in-place), v == 0, x == carry1
        // If there is an extra carry, propagate upwards
        mp_size_t pos = i + n + 1;
        while (c) {
            if (pos >= 2*n + 1) break;
            unsigned __int128 sum = (unsigned __int128)t[pos] + c;
            t[pos] = (mp_limb_t)sum;
            c = (mp_limb_t)(sum >> (sizeof(mp_limb_t)*8));
            ++pos;
        }
    }

    // At this point, t >> (n limbs) is the candidate result
    // Copy t[n .. 2n-1] to res
    for (mp_size_t i = 0; i < n; ++i) res[i] = t[i + n];

    // If res >= m, subtract once
    int ge = (mpn_cmp(res, m, n) >= 0);
    if (ge) {
        mpn_sub_n(res, res, m, n);
    }

    free(t);
    return 0;
}

int main(void) {
    // Example: 3000-bit numbers
    const unsigned int BITS = 3000;
    mpz_t A, B, M, R, T, Res_ref; // mpz high-level vars for testing
    mpz_inits(A, B, M, R, T, Res_ref, NULL);

    gmp_randstate_t st;
    gmp_randinit_default(st);

    mpz_urandomb(A, st, BITS);
    mpz_urandomb(B, st, BITS);

    mpz_urandomb(M, st, BITS);
    mpz_setbit(M, BITS-1);
    if (mpz_even_p(M)) mpz_add_ui(M, M, 1);

    // Convert to mpn limbs (n limbs)
    mp_size_t n_limbs = (BITS + GMP_NUMB_BITS - 1) / GMP_NUMB_BITS;

    mp_limb_t *a_limbs = calloc(n_limbs, sizeof(mp_limb_t));
    mp_limb_t *b_limbs = calloc(n_limbs, sizeof(mp_limb_t));
    mp_limb_t *m_limbs = calloc(n_limbs, sizeof(mp_limb_t));
    mp_limb_t *res_limbs = calloc(n_limbs, sizeof(mp_limb_t));

    size_t count;
    mpz_export(a_limbs, &count, -1, sizeof(mp_limb_t), 0, 0, A);
    // mpz_export writes least-significant limb first (little-endian) by default with -1
    // ensure count <= n_limbs
    mpz_export(b_limbs, &count, -1, sizeof(mp_limb_t), 0, 0, B);
    mpz_export(m_limbs, &count, -1, sizeof(mp_limb_t), 0, 0, M);

    // compute m0 inverse mod base
    mp_limb_t m0 = m_limbs[0];
    if ((m0 & 1) == 0) {
        fprintf(stderr, "Modulus must be odd for Montgomery\n");
        return 1;
    }
    mp_limb_t inv = inv_mod_base(m0);
    // We want m0inv_neg = (-inv) mod base such that m0 * m0inv_neg == -1 mod base
    mp_limb_t m0inv_neg = (mp_limb_t)(0 - inv);

    // Call montgomery multiplication (res = mont(a,b) in Montgomery domain)
    if (montgomery_mul_limbs(res_limbs, a_limbs, b_limbs, m_limbs, n_limbs, m0inv_neg) != 0) {
        fprintf(stderr, "malloc failed\n");
        return 1;
    }

    // Convert res_limbs back to mpz for verification
    mpz_t Res_low;
    mpz_init(Res_low);
    mpz_import(Res_low, n_limbs, -1, sizeof(mp_limb_t), 0, 0, res_limbs);

    // Verify against mpz arithmetic: compute (A*B*R^{-1}) mod M where R = base^n_limbs
    // Compute R = 2^{n_limbs * GMP_NUMB_BITS}
    mpz_set_ui(R, 0);
    mpz_setbit(R, n_limbs * GMP_NUMB_BITS);

    // canonical Montgomery product: Mont(a,b) = a*b*R^{-1} mod M
    // compute reference using mpz
    mpz_mul(T, A, B);
    mpz_mod(T, T, M);
    mpz_invert(Res_ref, R, M); // Res_ref = R^{-1} mod M (note: mpz_invert returns 0 if not invertible)
    if (mpz_cmp_ui(Res_ref, 0) == 0) {
        fprintf(stderr, "R not invertible mod M (unlikely)\n");
        return 1;
    }
    mpz_mul(T, T, Res_ref);
    mpz_mod(T, T, M);

    // compare Res_low and T
    if (mpz_cmp(Res_low, T) == 0) {
        printf("Montgomery multiplication (limb-level) matches reference.\n");
    } else {
        printf("Mismatch!\n");
        gmp_printf("res_limbs as mpz: %Zd\n", Res_low);
        gmp_printf("reference   : %Zd\n", T);
    }

    // cleanup
    free(a_limbs); free(b_limbs); free(m_limbs); free(res_limbs);
    mpz_clears(A, B, M, R, T, Res_ref, Res_low, NULL);
    gmp_randclear(st);
    return 0;
}
