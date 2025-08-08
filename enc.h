#pragma once

#include "shared.h"
#include "vfhe.h"

class Encryptor {
private:
i128 v, d, N, q;
vector_i128 sk;

public:
    Encryptor(i128 v_, i128 d_, i128 N_, i128 q_, vector_i128 sk_)
        : v(v_), d(d_), N(N_), q(q_), sk(sk_) {}
    // Encrypts an RLWE ciphertext of message m
    // m: message polynomial (poly), N: degree, sk: secret key, q: modulus, dth_pows: unused here
    rlwe encrypt_rlwe(const poly& m) {
        poly noise = poly(N);
        auto noise_vec = sample_noise_polynomial(N);
        for (size_t i = 0; i < N; ++i) noise.set(i, noise_vec[i]);

        poly a = poly(N);
        auto a_vec = sample_random_polynomial(N, q);
        for (size_t i = 0; i < N; ++i) a.set(i, a_vec[i]);

        // poly_mult_mod(a, sk, q)
        poly ask(N);
        for (size_t i = 0; i < N; ++i) {
            i128 sum = 0;
            for (size_t j = 0; j < N; ++j) {
                sum += a.get(j) * sk[mod_((i - j + N), N)];
            }
            ask.set(i, mod_(sum, q));
        }

        poly b(N);
        for (size_t i = 0; i < N; ++i) {
            i128 val = m.get(i) + noise.get(i) + ask.get(i);
            b.set(i, mod_(val, q));
        }

        rlwe res(2, N);
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

    // Encrypts an RGSW ciphertext of message M (poly)
    rgsw encrypt_rgsw(const poly& M) {
        // Compute v powers: v^0, v^1, ..., v^{d-1}
        vector_i128 v_powers(d);
        for (i128 i = 0; i < d; ++i) {
            v_powers[i] = 1;
            for (i128 j = 0; j < i; ++j) {
                v_powers[i] = mod_((v_powers[i] * v), q);
            }
        }

        // Build G matrix: 2 x 2d, each row is [v_powers, 0...], [0..., v_powers]
        std::vector<vector_i128> G(2, vector_i128(2 * d, 0));
        for (i128 i = 0; i < d; ++i) {
            G[0][i] = v_powers[i];
            G[1][d + i] = v_powers[i];
        }

        // Encryptions of zero
        std::vector<rlwe> encs_of_zero(2 * d);
        poly zero_poly(N);
        for (i128 i = 0; i < 2 * d; ++i) {
            encs_of_zero[i] = encrypt_rlwe(zero_poly);
        }

        // Compute M * G and add to RLWE encryptions
        // NOTE needs modification for
        // FIXME
        for (i128 i = 0; i < N; ++i) {
            for (i128 j = 0; j < 2 * d; ++j) {
                // Add M[i] * G[row][j] to b poly of RLWE
                for (i128 k = 0; k < 2; k++) {
                    i128 val = mod_(M.get(i) * G.at(k).at(j), q);
                    rlwe& ct = encs_of_zero.at(j);
                    poly& rlwe_component = ct.get_poly(k);
                    rlwe_component.set(i, mod_(rlwe_component.get(i) + val, q));
                }
            }
        }

        // Pack RLWE encryptions into RGSW
        rgsw res(2 * d, 2, N);
        for (i128 i = 0; i < 2 * d; ++i) {
            res.set(i, encs_of_zero[i]);
        }
        return res;
    }
    // takes an array2d<i128> and returns an encrypted rgsw_mat
    rgsw_mat encrypt_rgsw_mat(const array2d<i128>& mat) {
        size_t rows = mat.n_rows();
        size_t cols = mat.n_cols();
        rgsw_mat res(rows, cols, 2 * d, 2, N);
        for (size_t i = 0; i < rows; ++i) {
            for (size_t j = 0; j < cols; ++j) {
                poly p(N);
                for (size_t k = 0; k < N; ++k) {
                    p.set(k, mat.get(i, j)); // fill all coeffs with mat(i, j)
                }
                rgsw enc_rgsw = encrypt_rgsw(p);
                res.set(i, j, enc_rgsw);
            }
        }
        return res;
    }

    // Encodes an RGSW plaintext (not encryption, just encoding)
    rgsw encode_rgsw(const poly& M) {
        // Compute v powers: v^0, v^1, ..., v^{d-1}
        vector_i128 v_powers(d);
        for (i128 i = 0; i < d; ++i) {
            v_powers[i] = 1;
            for (i128 j = 0; j < i; ++j) {
                v_powers[i] = mod_(v_powers[i] * v, q);
            }
        }

        // Build G matrix: 2 x 2d, each row is [v_powers, 0...], [0..., v_powers]
        std::vector<vector_i128> G(2, vector_i128(2 * d, 0));
        for (i128 i = 0; i < d; ++i) {
            G[0][i] = v_powers[i];
            G[1][d + i] = v_powers[i];
        }

        // Create rgsw object
        rgsw res(2 * d, 2, N);
        for (i128 j = 0; j < 2 * d; ++j) {
            rlwe ct(2, N);
            for (i128 row = 0; row < 2; ++row) {
                poly p(N);
                for (i128 i = 0; i < N; ++i) {
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
