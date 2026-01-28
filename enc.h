#pragma once

#include "shared.h"
#include "vfhe.h"
#include "gmpxx.h"

poly sample_discrete_gaussian(size_t N, double mu = 3.2, double sigma = 19.2) {
    poly result(N);
    std::random_device rd;
    std::mt19937 gen(rd());
    std::normal_distribution<double> dist(mu, sigma);
    for (size_t i = 0; i < N; i++) {
        bigz res = static_cast<long int>(std::round(dist(gen)));
        if (res < 0)
            res += FIELD_MODULUS;
        result[i] = static_cast<bigz>(res);
    }
    #ifdef DEBUG1_ON
        for (size_t i = 0; i < N; i++)
            result[i] = (i % (N - 1) == 0) ? 0 : 0;
    #endif
    return result;
}

poly sample_secret_key(size_t N) {
    // Sample a secret key for the RGSW scheme.
    // Each entry is -1, 0, or 1, with probabilities 0.25, 0.5, 0.25 respectively.
    poly s(N);
    std::random_device rd;
    std::mt19937 gen(rd());
    std::discrete_distribution<int> dist({0.25, 0.5, 0.25});
    for (size_t i = 0; i < N; i++) {
        int val = dist(gen);
        // TODO uncomment and deal with potential check_val_bounds errors
        // if (val == 0) s[i] = -1;
        if (val == 0) s[i] = 1;
        else if (val == 1) s[i] = 0;
        // else s[i] = FIELD_MODULUS - 1; // FIXME
        else s[i] = 1;
    }
    #ifdef DEBUG1_ON
    for (size_t i = 0; i < N; i++)
        s[i] = 1;
    #endif
    return s;
}
// Sample a random polynomial of degree N-1 with coefficients in the range [0, q).
poly sample_random_polynomial(size_t N, const bigz& q) {
    poly ply(N);
    // std::random_device rd;
    // std::mt19937 gen(rd());
    // std::uniform_int_distribution<bigz> dist(0, q - 1);
    bigz rop;
    gmp_randstate_t state;
    gmp_randinit_default(state);
    bigz seed = 42;
    gmp_randseed(state, seed.get_mpz_t());
    // mp_bitcnt_t n = log2_mpz(q) + 1;
    for (size_t i = 0; i < N; i++) {
        // poly[i] = dist(gen);
        // mpz_urandomb(rop, state, n);
        mpz_urandomm(rop.get_mpz_t(), state, q.get_mpz_t());
        ply[i] = bigz(rop);
    }
    gmp_randclear(state);
    #ifdef DEBUG1_ON
        for (size_t i = 0; i < N; i++)
            ply[i] = 1;
    #endif
    return ply;
}

// Sample a noise polynomial of degree N-1 with coefficients from the discrete Gaussian distribution.
// poly sample_noise_polynomial(size_t N, double mu = 3.2, double sigma = 19.2) {
poly sample_noise_polynomial(size_t N) {
    // return sample_discrete_gaussian(N, mu, sigma);
    return sample_secret_key(N); // FIXME
}

class Encryptor {
private:
    bigz v;
    size_t d;
    size_t N;
    bigz q;
    poly sk;
public:
    Encryptor(bigz v_, size_t d_, size_t N_, bigz q_)
        : v(v_), d(d_), N(N_), q(q_), sk{sample_secret_key(N).to_eval_form(N, poly::k - 1)} { }
    // Encrypts an RLWE ciphertext of message m
    // m: message polynomial (poly), N: degree, sk: secret key, q: modulus, dth_pows: unused here
    rlwe encrypt_rlwe(const poly& m) {
        poly noise = sample_noise_polynomial(N);
        poly a = sample_random_polynomial(N, q);
        poly b = m + noise + a * sk;
        rlwe res(N_POLYS_IN_RLWE);
        res.set(0, b);
        res.set(1, a);
        return res;
    }

    poly decrypt_rlwe(const rlwe& ctx) const {
        ASSERT(ctx.n_polys() == 2);
        const poly& b = ctx.get_poly(0);
        const poly& a = ctx.get_poly(1);
        poly m = b - a * sk;
        return m;
    }

    rlwe_vec encrypt_rlwe_vec(const vector_bigz& vec) {
        rlwe_vec res(vec.size());
        // #pragma omp parallel for schedule(static) num_threads(N_THREADS) // FIXME
        for (size_t i = 0; i < vec.size(); i++) {
            poly p{N};
            p.set(0, vec.at(i)); // set the first coefficient to the value
            rlwe r = encrypt_rlwe(p);
            res.set(i, r);
        }
        return res;
    }

    vector_bigz decrypt_rlwe_vec(const rlwe_vec& ctxs) const {
        vector_bigz ptxs(ctxs.size());
        for (size_t i = 0; i < ctxs.size(); i++) {
            poly m = decrypt_rlwe(ctxs.get(i));
            ptxs.at(i) = m.get(0);
        }
        return ptxs;
    }
    rgsw encrypt_rgsw(const poly& M) {
        // Compute v powers: v^0, v^1, ..., v^{d-1}
        // TODO move into constructor
        TIMING(static auto start = std::chrono::high_resolution_clock::now();)
        vector_bigz v_powers(d);
        v_powers[0] = 1;
        for (size_t i = 1; i < d; i++)
            v_powers[i] = v_powers[i - 1] * v % FIELD_MODULUS;
        TIMING(static auto end = std::chrono::high_resolution_clock::now();)
        TIMING(timing.v_pows += end - start;)

        // Build (2x2d) G matrix: with rows: [v_powers, 0...], [0..., v_powers]
        TIMING(start = std::chrono::high_resolution_clock::now();)
        rgsw G(2 * d, N_POLYS_IN_RLWE, N);
        for (size_t i = 0; i < d; i++) {
            G.get_rlwe(i).get_poly(0).set(0, mod_(M.get(0) * v_powers[i], FIELD_MODULUS));
            G.get_rlwe(i + d).get_poly(1).set(0, mod_(M.get(0) * v_powers[i], FIELD_MODULUS));
        }
        TIMING(end = std::chrono::high_resolution_clock::now();)
        TIMING(timing.g_mat += end - start;)

        // Encryptions of zero
        TIMING(start = std::chrono::high_resolution_clock::now();)
        rgsw encs_of_zero(2 * d);
        poly zero_poly(N);
        // #pragma omp parallel for schedule(static) num_threads(N_THREADS) // FIXME
        for (size_t i = 0; i < 2 * d; i++)
            encs_of_zero.set(i, encrypt_rlwe(zero_poly));
        TIMING(end = std::chrono::high_resolution_clock::now();)
        TIMING(timing.enc_zero += end - start;)

        // rgsw res = G + encs_of_zero;
        TIMING(start = std::chrono::high_resolution_clock::now();)
        G += encs_of_zero;
        TIMING(end = std::chrono::high_resolution_clock::now();)
        TIMING(timing.enc_rgsw_add += end - start;)
        return G;
    }
    // takes an array2d<bigz> and returns an encrypted rgsw_mat
    rgsw_mat encrypt_rgsw_mat(const array2d<bigz>& mat) {
        size_t rows = mat.n_rows();
        size_t cols = mat.n_cols();
        rgsw_mat res(rows, cols); //, 2 * d, N_POLYS_IN_RLWE, N);
        for (size_t i = 0; i < rows; i++) {
            for (size_t j = 0; j < cols; j++) {
                poly p(N);
                // Set only the poly's constant coeff
                // bigz val = mat.get(i, j);
                p.set(0, mod_(mat.get(i, j), FIELD_MODULUS));
                res.set(i, j, encrypt_rgsw(p));
            }
        }
        return res;
    }

    // Encodes an RGSW plaintext (not encryption, just encoding)
    flat_rgsw encode_flat_rgsw(const poly& M0, const poly& M1) const {
        // Compute v powers: v^0, v^1, ..., v^{d-1}
        // TODO move into Encryptor constructor
        vector_bigz v_powers(d);
        v_powers[0] = 1;
        for (size_t i = 1; i < d; i++) {
            v_powers.at(i) = v * v_powers.at(i - 1);
            v_powers.at(i) = mod_(v_powers.at(i), q);
        }

        // Create flat_rgsw object
        flat_rgsw res(2 * d);
        for (size_t j = 0; j < d; j++) {
            res.set(j, M0 * v_powers.at(j));
            res.set(j + d, M1 * v_powers.at(j));
        }
        return res;
    }
};
