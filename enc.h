#pragma once

#include "shared.h"
#include "vfhe.h"

vector_i128 sample_discrete_gaussian(size_t N, double mu = 3.2, double sigma = 19.2) {
    vector_i128 result(N);
    std::random_device rd;
    std::mt19937 gen(rd());
    std::normal_distribution<double> dist(mu, sigma);
    for (size_t i = 0; i < N; i++) {
        i128 res = static_cast<i128>(std::round(dist(gen)));
        result[i] = mod_(res, FIELD_MODULUS);
    }
    #ifdef DEBUG1_ON
        for (size_t i = 0; i < N; i++)
            result.at(i) = 1;
    #endif
    return result;
}

vector_i128 sample_secret_key(size_t N) {
    // Sample a secret key for the RGSW scheme.
    // Each entry is -1, 0, or 1, with probabilities 0.25, 0.5, 0.25 respectively.
    vector_i128 s(N);
    std::random_device rd;
    std::mt19937 gen(rd());
    std::discrete_distribution<int> dist({0.25, 0.5, 0.25});
    for (size_t i = 0; i < N; i++) {
        int val = dist(gen);
        // TODO uncomment and deal with potential check_val_bounds errors
        // if (val == 0) s[i] = -1;
        if (val == 0) s[i] = 1;
        else if (val == 1) s[i] = 0;
        else s[i] = 1;
    }
    #ifdef DEBUG1_ON
    for (size_t i = 0; i < N; i++)
        s.at(i) = 1;
    #endif
    return s;
}
// Sample a random polynomial of degree N-1 with coefficients in the range [0, q).
vector_i128 sample_random_polynomial(size_t N, i128 q) {
    vector_i128 poly(N);
    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_int_distribution<i128> dist(0, q - 1);
    for (size_t i = 0; i < N; i++) {
        poly[i] = dist(gen);
    }
    #ifdef DEBUG1_ON
        for (size_t i = 0; i < N; i++)
            poly.at(i) = 1;
    #endif
    return poly;
}

// Sample a noise polynomial of degree N-1 with coefficients from the discrete Gaussian distribution.
vector_i128 sample_noise_polynomial(size_t N, double mu = 3.2, double sigma = 19.2) {
    return sample_discrete_gaussian(N, mu, sigma);
}

class Encryptor {
private:
    u32 v, d;
    size_t N;
    i128 q;

public:
    poly sk;
    Encryptor(u32 v_, u32 d_, size_t N_, i128 q_)
        : v(v_), d(d_), N(N_), q(q_), sk(N) {
            vector_i128 sk_ = sample_secret_key(N);
            for (size_t i = 0; i < N; i++)
                sk.set(i, sk_.at(i));
        }
    // Encrypts an RLWE ciphertext of message m
    // m: message polynomial (poly), N: degree, sk: secret key, q: modulus, dth_pows: unused here
    rlwe encrypt_rlwe(const poly& m) {
        poly noise = poly(N);
        auto noise_vec = sample_noise_polynomial(N);
        for (size_t i = 0; i < N; i++)
            // noise.set(i, noise_vec[i], true);
            noise.set(i, noise_vec[i]);
        // for (size_t i = 0; i < 10; i++) {
        //     auto val = noise.get(i);
        //     (i < 9) ? (std::cout << i128str(val) << ", ") : (std::cout << i128str(val) << "\n");
        // }

        poly a = poly(N);
        auto a_vec = sample_random_polynomial(N, q);
        for (size_t i = 0; i < N; i++)
            a.set(i, a_vec[i]);

        poly ask = a * sk;

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
        // Compute the plaintext polynomial
        poly m = b - a * sk;
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
        v_powers[0] = 1;
        for (size_t i = 1; i < d; i++) {
            v_powers[i] = v_powers[i - 1] * v % FIELD_MODULUS;
        }

        // Build G matrix: 2d x 2, each col is [v_powers, 0...], [0..., v_powers]
        std::vector<vector_i128> G(2 * d, vector_i128(2, 0));
        for (size_t i = 0; i < d; i++) {
            G[i][0] = v_powers[i];
            G[d + i][1] = v_powers[i];
        }

        // Encryptions of zero
        std::vector<rlwe> encs_of_zero(2 * d);
        poly zero_poly(N);
        for (size_t i = 0; i < 2 * d; i++) {
            encs_of_zero[i] = encrypt_rlwe(zero_poly);
        }

        // Compute M * G and add to RLWE encryptions
        for (size_t i = 0; i < 2 * d; i++) {
            // Add M[i] * G[row][j] to b poly of RLWE
            std::cout << "rlwe[" << i << "]:";
            for (size_t j = 0; j < N_POLYS_IN_RLWE; j++) {
                poly val = M * G.at(i).at(j);
                if (j == 0)
                    std::cout << "(" << i128str(val.get(0)) << ", ";
                else
                    std::cout << i128str(val.get(0)) << ")";
                rlwe& ct = encs_of_zero.at(i);
                ct.set(j, ct.get(j) + val);
            }
            std::cout << "\n";
        }

        // Pack RLWE encryptions into RGSW
        rgsw res(2 * d, N_POLYS_IN_RLWE, N);
        for (size_t i = 0; i < 2 * d; i++) {
            res.set(i, encs_of_zero[i]);
        }

        for (size_t i = 0; i < 2 * d; i++) {
            // Add M[i] * G[row][j] to b poly of RLWE
            std::cout << "rgsw[" << i << "]:";
            for (size_t j = 0; j < N_POLYS_IN_RLWE; j++) {
                poly val = res.get(i).get(j);
                if (j == 0)
                    std::cout << "(" << i128str(val.get(0)) << ", ";
                else
                    std::cout << i128str(val.get(0)) << ")";
            }
            std::cout << "\n";
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
