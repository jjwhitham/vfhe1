#pragma once

#include "shared.h"
#include "vfhe.h"

// // Find prime q: q = 1 + 2*N*k, k >= 1, q % (2*N) == 1, q is prime
// i128 find_q(i128 N, i128 q_init) {
//     int i = 0;
//     i128 k = ((q_init - 1) / (2 * N)) + 1;
//     i128 q = 1 + 2 * N * k;
//     while (i < 100) {
//         if ((q % (2 * N)) == 1 && is_prime(q)) {
//             return q;
//         }
//         i++;
//         k++;
//         q = 1 + 2 * N * k;
//     }
//     // If not found, return -1 or throw
//     throw std::runtime_error("No suitable prime q found.");
// }

// // Generation of a cyclic group of prime order
// // Find prime p: p = k * q + 1, k >= 2
// // Given a prime q, find another prime p such that q | (p-1)
// std::pair<i128, i128> get_parent_group(i128 q) {
//     if (!is_prime(q)) {
//         throw std::invalid_argument("q must be a prime number.");
//     }
//     i128 k = 2;
//     while (!is_prime(k * q + 1)) {
//         k++;
//         if (k == 1000000) {
//             throw std::runtime_error("No suitable prime found for the given order.");
//         }
//     }
//     i128 p = k * q + 1;
//     return {p, k};
// }

// // The cyclic group Z/pZ* (of order p - 1) contains a cyclic subgroup of order q.
// // Find a generator of this cyclic subgroup.
// i128 get_generator_of_prime_order_group(i128 p, i128 q, i128 k) {
//     // Randomly pick a number from [2, q-1] and check if it is a generator of the cyclic subgroup of order q.
//     for (i128 h = 2; h < q; ++h) {
//         i128 g = mod_pow(h, k, p);
//         if (g != 1) {
//             // Candidate is a generator of order q
//             if (mod_pow(g, q, p) == 1) {
//                 return g;
//             }
//         }
//     }
//     throw std::runtime_error("No generator found for the given prime order.");
// }

// // For a given finite field of size q (prime), generate a cyclic group of prime order q.
// // q is a prime number ~= 2**54
// // Example: q = (5 * 7 * 62829235873) * 2^13 + 1 (form: m * 2^k + 1 for NTT)
// // q = 18014398509506561
// // q = 17
// // q = 1073750017
// std::tuple<i128, i128, i128> generate_field_and_group_params() {
//     i128 q = find_q(1 << 12, 1 << 12);
//     auto [p, k] = get_parent_group(q);
//     i128 g = get_generator_of_prime_order_group(p, q, k);
//     std::cout << "p: " << p << std::endl;
//     std::cout << "q: " << q << std::endl;
//     std::cout << "k: " << k << std::endl;
//     std::cout << "g: " << g << std::endl << std::endl;
//     // std::cout << "Generated field with prime " << q << " and cyclic group of order " << p-1 << " with generator " << g << "." << std::endl;
//     return {g, q, p};
// }

class Encryptor {
private:
    u32 v, d;
    size_t N;
    i128 q;
    vector_i128 sk;

public:
    Encryptor(u32 v_, u32 d_, size_t N_, i128 q_, vector_i128 sk_)
        : v(v_), d(d_), N(N_), q(q_), sk(sk_) {}
    // Encrypts an RLWE ciphertext of message m
    // m: message polynomial (poly), N: degree, sk: secret key, q: modulus, dth_pows: unused here
    rlwe encrypt_rlwe(const poly& m) {
        poly noise = poly(N);
        auto noise_vec = sample_noise_polynomial(N);
        for (size_t i = 0; i < N; i++)
            // noise.set(i, noise_vec[i], true);
            noise.set(i, noise_vec[i]);

        poly a = poly(N);
        auto a_vec = sample_random_polynomial(N, q);
        for (size_t i = 0; i < N; i++)
            a.set(i, a_vec[i]);

        // make poly from sk
        // FIXME polys all the way through
        poly sk_poly(N);
        for (size_t i = 0; i < N; i++)
            // sk_poly.set(i, sk.at(i), true);
            sk_poly.set(i, sk.at(i));
        poly ask = a * sk_poly;

        poly b(N);
        for (size_t i = 0; i < N; i++) {
            i128 val = m.get(i) + noise.get(i) + ask.get(i);
            b.set(i, mod_(val, q));
        }

        rlwe res(N_POLYS_IN_RLWE, N);
        res.set(0, b);
        res.set(1, a);
        return res;
    }

    rlwe_vec encrypt_rlwe_vec(const vector_i128& vec) {
        rlwe_vec res(vec.size());
        for (size_t i = 0; i < vec.size(); i++) {
            poly p(N);
            p.set(0, vec.at(i)); // set the first coefficient to the value
            rlwe r = encrypt_rlwe(p);
            res.set(i, r);
        }
        return res;
    }

    poly decrypt_rlwe(const rlwe& ctx) const {
        ASSERT(ctx.n_polys() == 2);
        poly b = ctx.get_poly(0);
        poly a = ctx.get_poly(1);

        // Compute the inverse of the secret key
        poly sk_poly(N);
        for (size_t i = 0; i < N; i++)
            sk_poly.set(i, sk.at(i));

        // Compute the plaintext polynomial
        poly m = b - a * sk_poly;
        return m;
    }

    vector_i128 decrypt_rlwe_vec(const rlwe_vec& ctxs) const {
        vector_i128 ptxs(ctxs.size());
        for (size_t i = 0; i < ctxs.size(); i++) {
            poly m = decrypt_rlwe(ctxs.get(i));
            ptxs.at(i) = m.get(0);
        }
        return ptxs;
    }
    // Encrypts an RGSW ciphertext of message M (poly)
    rgsw encrypt_rgsw(const poly& M) {
        // Compute v powers: v^0, v^1, ..., v^{d-1}
        vector_i128 v_powers(d);
        for (size_t i = 0; i < d; i++) {
            v_powers[i] = 1;
            for (size_t j = 0; j < i; j++) {
                v_powers[i] = mod_((v_powers[i] * v), q);
            }
        }

        // Build G matrix: 2 x 2d, each row is [v_powers, 0...], [0..., v_powers]
        std::vector<vector_i128> G(2, vector_i128(2 * d, 0));
        for (size_t i = 0; i < d; i++) {
            G[0][i] = v_powers[i];
            G[1][d + i] = v_powers[i];
        }

        // Encryptions of zero
        std::vector<rlwe> encs_of_zero(2 * d);
        // FIXME is this correct? I'm using the same zero poly for each call
        poly zero_poly(N);
        for (size_t i = 0; i < 2 * d; i++) {
            encs_of_zero[i] = encrypt_rlwe(zero_poly);
        }

        // Compute M * G and add to RLWE encryptions
        // NOTE needs modification for
        // FIXME
        for (size_t i = 0; i < N; i++) {
            for (size_t j = 0; j < 2 * d; j++) {
                // Add M[i] * G[row][j] to b poly of RLWE
                for (size_t k = 0; k < 2; k++) {
                    i128 val = mod_(M.get(i) * G.at(k).at(j), q);
                    rlwe& ct = encs_of_zero.at(j);
                    poly& rlwe_component = ct.get_poly(k);
                    rlwe_component.set(i, mod_(rlwe_component.get(i) + val, q));
                }
            }
        }

        // Pack RLWE encryptions into RGSW
        rgsw res(2 * d, N_POLYS_IN_RLWE, N);
        for (size_t i = 0; i < 2 * d; i++) {
            res.set(i, encs_of_zero[i]);
        }
        return res;
    }
    // takes an array2d<i128> and returns an encrypted rgsw_mat
    rgsw_mat encrypt_rgsw_mat(const array2d<i128>& mat) {
        size_t rows = mat.n_rows();
        size_t cols = mat.n_cols();
        rgsw_mat res(rows, cols, 2 * d, N_POLYS_IN_RLWE, N);
        for (size_t i = 0; i < rows; i++) {
            for (size_t j = 0; j < cols; j++) {
                poly p(N);
                // Set only the poly's constant coeff
                p.set(0, mod_(mat.get(i, j), FIELD_MODULUS));
                rgsw enc_rgsw = encrypt_rgsw(p);
                res.set(i, j, enc_rgsw);
            }
        }
        return res;
    }

    // Encodes an RGSW plaintext (not encryption, just encoding)
    rgsw encode_rgsw(const poly& M) const {
        // Compute v powers: v^0, v^1, ..., v^{d-1}
        vector_i128 v_powers(d);
        for (size_t i = 0; i < d; i++) {
            v_powers[i] = 1;
            for (size_t j = 0; j < i; j++) {
                v_powers[i] = mod_(v_powers[i] * v, q);
            }
        }

        // Build G matrix: 2 x 2d, each row is [v_powers, 0...], [0..., v_powers]
        std::vector<vector_i128> G(2, vector_i128(2 * d, 0));
        for (size_t i = 0; i < d; i++) {
            G[0][i] = v_powers[i];
            G[1][d + i] = v_powers[i];
        }

        // Create rgsw object
        rgsw res(2 * d, 2, N);
        for (size_t j = 0; j < 2 * d; j++) {
            rlwe ct(2, N);
            for (size_t row = 0; row < 2; ++row) {
                poly p(N);
                for (size_t i = 0; i < N; i++) {
                    i128 val = mod_(M.get(i) * G[row][j], q);
                    p.set(i, val);
                }
                ct.set(row, p);
            }
            res.set(j, ct);
        }
        return res;
    }
};
