// TODO write more tests
// TODO rgsw code
    // TODO PRNG
    // TODO decomposition
    // TODO enc/dec
// TODO control code
// TODO verification code
// TODO proof code (linear/dynamic checks)
// TODO NTT code
// TODO homomorphic hash variant
// TODO parallellise
// TODO cyclic group code
// TODO replace *this and (*this).member_func() with this and this->member_func() ???
/*
Setup:
rF, rG, sH
- veri_vec * rgsw_vec
    - scalar * rgsw

Verfication:
rx' = rFx + rGx
su = sHu
*/


#include "shared.h"
#include "enc.h"
#include "vfhe.h"
#include <ranges>


vector_i128 eval_poly_pows(size_t n, i128 base, i128 q) {
    vector_i128 res(n);
    res.at(0) = 1; // base^0 = 1
    for (size_t i = 1; i < n; i++) {
        res.at(i) = mod_(res.at(i - 1) * base, q);
    }
    return res;
}

std::tuple<eval_key, veri_key> compute_eval_and_veri_keys(
    rgsw_mat F_ctx, rgsw_mat G_bar_ctx, rgsw_mat R_bar_ctx, rgsw_mat H_bar_ctx,
    vector_i128 r_0, vector_i128 r_1, vector_i128 s,
    i128 rho_0, i128 rho_1, i128 alpha_0, i128 alpha_1, i128 gamma_0, i128 gamma_1,
    i128 g, i128 d, i128 q, i128 p, i128 N, vector_i128 eval_pows, Encryptor enc
    // i128 g, i128 d, i128 q, i128 p, i128 v, i128 N, vector_i128 eval_pows, Encryptor enc
) {
    // size_t n = F_ctx.n_rows();
    // size_t m = H_bar_ctx.n_rows();

    // Anonymous function for cyclic exponentiation: (g^e1)^e2 mod p
    auto cyclic_exp = [](i128 base, i128 exp, i128 q, i128 p) -> i128 {
        // Compute base^exp mod p
        i128 result = 1;
        base = mod_(base, p);
        exp = mod_(exp, q);
        while (exp > 0) {
            if (exp & 1)
                result = mod_(result * base, p);
            exp >>= 1;
            base = mod_(base * base, p);
        }
        return result;
    };

    // Converts a vector v to a vector where each element is g^v[i] mod p
    auto convert_vec_to_cyclic = [&](i128 g, const vector_i128& v, i128 q, i128 p) -> vector_i128 {
        vector_i128 res(v.size());
        for (size_t i = 0; i < v.size(); ++i) {
            res[i] = cyclic_exp(g, v[i], q, p);
        }
        return res;
    };

    vector_i128 gr_0 = convert_vec_to_cyclic(g, r_0, q, p);
    vector_i128 gr_1 = convert_vec_to_cyclic(g, r_1, q, p);
    vector_i128 gr_rho_0 = convert_vec_to_cyclic(g, scalar_vec_mult(rho_0, r_0, q), q, p); // Example for first element
    vector_i128 gr_rho_1 = convert_vec_to_cyclic(g, scalar_vec_mult(rho_1, r_1, q), q, p);

    // Anonymous function for vector-matrix multiplication: vector_i128 * rgsw_mat
    auto vec_mat_mult = [](const vector_i128& vec, const rgsw_mat& mat) -> rgsw_vec {
        size_t rows = mat.n_rows();
        size_t cols = mat.n_cols();
        assert(vec.size() == rows);
        rgsw_vec res(cols, mat.get_n_rlwes(), mat.get_n_polys(), mat.get_n_coeffs());
        for (size_t j = 0; j < cols; ++j) {
            rgsw sum(mat.get_n_rlwes(), mat.get_n_polys(), mat.get_n_coeffs());
            for (size_t i = 0; i < rows; ++i) {
                sum = sum + (mat.get(i, j) * vec[i]);
            }
            res.set(j, sum);
        }
        return res;
    };

    // Hashes a poly using eval_pows and q
    auto hash_poly = [](const poly& p, const vector_i128& eval_pows, i128 q) -> i128 {
        i128 sum = 0;
        for (size_t i = 0; i < p.size() && i < eval_pows.size(); ++i) {
            sum = mod_(sum + mod_(p.get(i) * eval_pows[i], q), q);
        }
        return sum;
    };

    // Hashes an rgsw_vec using eval_pows and q
    auto hash_rgsw_vector = [&](const rgsw_vec& rgsw_v, const vector_i128& eval_pows, i128 q) -> rgsw_vec {
        rgsw_vec hashed(rgsw_v.size(), rgsw_v.get_n_rlwes(), rgsw_v.get_n_polys(), rgsw_v.get_n_coeffs());
        for (size_t i = 0; i < rgsw_v.size(); ++i) {
            const rgsw& rgsw_elem = rgsw_v.get(i);
            rgsw hashed_rgsw(rgsw_elem.size(), rgsw_elem.n_polys(), rgsw_elem.n_coeffs());
            for (size_t j = 0; j < rgsw_elem.size(); ++j) {
                const rlwe& rlwe_elem = rgsw_elem.get(j);
                rlwe hashed_rlwe(rlwe_elem.size(), rlwe_elem.get(0).size());
                for (size_t k = 0; k < rlwe_elem.size(); ++k) {
                    const poly& poly_elem = rlwe_elem.get(k);
                    i128 hash_val = hash_poly(poly_elem, eval_pows, q);
                    poly hashed_poly(poly_elem.size());
                    for (size_t l = 0; l < poly_elem.size(); ++l) {
                        hashed_poly.set(l, hash_val);
                    }
                    hashed_rlwe.set(k, hashed_poly);
                }
                hashed_rgsw.set(j, hashed_rlwe);
            }
            hashed.set(i, hashed_rgsw);
        }
        return hashed;
    };

    rgsw_vec rF_0 = vec_mat_mult(r_0, F_ctx);
    rgsw_vec rF_1 = vec_mat_mult(r_1, F_ctx);
    rF_0 = hash_rgsw_vector(rF_0, eval_pows, q);
    rF_1 = hash_rgsw_vector(rF_1, eval_pows, q);

    // Make rgsw_vec for r_0 and r_1
    rgsw_vec r_0_rgsw(r_0.size(), 2 * d, 2, N), r_1_rgsw(r_1.size(), 2 * d, 2, N);
    for (size_t i = 0; i < r_0.size(); ++i) {
        poly p(N);
        p.set(0, r_0[i]);
        r_0_rgsw.set(i, enc.encode_rgsw(p));
    }
    for (size_t i = 0; i < r_1.size(); ++i) {
        poly p(N);
        p.set(0, r_1[i]);
        r_1_rgsw.set(i, enc.encode_rgsw(p));
    }
    assert(r_0_rgsw.size() == rF_0.size());
    // assert(r_0_rgsw[0].rows() == rF_0[0].rows()); // FIXME

    rgsw_vec rF_0_r_1(rF_0.size()), rF_1_r_0(rF_1.size());
    for (size_t i = 0; i < rF_0.size(); ++i)
        rF_0_r_1.set(i, rF_0.get(i) - r_1_rgsw.get(i));
    for (size_t i = 0; i < rF_1.size(); ++i)
        rF_1_r_0.set(i, rF_1.get(i) - r_0_rgsw.get(i));

    // = convert_vec_to_cyclic_enc(g, rF_0_r_1, d, q, p);
    rgsw_vec grFr_0 = rF_0_r_1.pow();
    // convert_vec_to_cyclic_enc(g, rF_1_r_0, d, q, p);
    rgsw_vec grFr_1 = rF_1_r_0.pow();

    // convert_vec_to_cyclic_enc(g, rF_0_r_1 * alpha_1, d, q, p);
    rgsw_vec grFr_alpha_0 = (rF_0_r_1 * alpha_1).pow();
    // convert_vec_to_cyclic_enc(g, rF_1_r_0 * alpha_0, d, q, p);
    rgsw_vec grFr_alpha_1 = (rF_1_r_0 * alpha_0).pow();

    rgsw_vec rG_0 = vec_mat_mult(r_0, G_bar_ctx);
    rgsw_vec rG_1 = vec_mat_mult(r_1, G_bar_ctx);
    rgsw_vec rR_0 = vec_mat_mult(r_0, R_bar_ctx);
    rgsw_vec rR_1 = vec_mat_mult(r_1, R_bar_ctx);
    rG_0 = hash_rgsw_vector(rG_0, eval_pows, q);
    rG_1 = hash_rgsw_vector(rG_1, eval_pows, q);
    rR_0 = hash_rgsw_vector(rR_0, eval_pows, q);
    rR_1 = hash_rgsw_vector(rR_1, eval_pows, q);

    rgsw_vec sH = vec_mat_mult(s, H_bar_ctx);
    sH = hash_rgsw_vector(sH, eval_pows, q);

    rgsw_vec sH_r_1(sH.size()), sH_r_0(sH.size());
    for (size_t i = 0; i < sH.size(); ++i)
        sH_r_1.set(i, sH.get(i) - r_1_rgsw.get(i));
    for (size_t i = 0; i < sH.size(); ++i)
        sH_r_0.set(i, sH.get(i) - r_0_rgsw.get(i));

    // rgsw_vec gsHr_0 = convert_vec_to_cyclic_enc(g, mod_rgsw_vec(sH_r_1, q), d, q, p);
    rgsw_vec gsHr_0 = sH_r_1.pow();
    // rgsw_vec gsHr_1 = convert_vec_to_cyclic_enc(g, mod_rgsw_vec(sH_r_0, q), d, q, p);
    rgsw_vec gsHr_1 = sH_r_0.pow();

    // rgsw_vec gsHr_gamma_0 = convert_vec_to_cyclic_enc(g, mod_rgsw_vec(sH_r_1 * gamma_1, q), d, q, p);
    rgsw_vec gsHr_gamma_0 = (sH_r_1 * gamma_1).pow();
    // rgsw_vec gsHr_gamma_1 = convert_vec_to_cyclic_enc(g, mod_rgsw_vec(sH_r_0 * gamma_0, q), d, q, p);
    rgsw_vec gsHr_gamma_1 = (sH_r_0 * gamma_0).pow();

    eval_key ek {
        gr_0,
        grFr_0,
        gsHr_0,
        gr_rho_0,
        grFr_alpha_0,
        gsHr_gamma_0,
        gr_1,
        grFr_1,
        gsHr_1,
        gr_rho_1,
        grFr_alpha_1,
        gsHr_gamma_1
    };

    veri_key vk {
        s,
        rG_0,
        rG_1,
        rR_0,
        rR_1,
        rho_0,
        rho_1,
        alpha_0,
        alpha_1,
        gamma_0,
        gamma_1
    };

    return std::make_tuple(ek, vk);
}

// TODO
// g, q, p = generate_field_and_group_params()

void run_control_loop() {
    Params pms;
    using matrix_i128 = array2d<i128>;

    // Helper function for scalar-vector multiplication (mod q)
    auto scalar_vec_mult = [](i128 scalar, const vector_i128& vec, i128 q) -> vector_i128 {
        vector_i128 result(vec.size());
        for (size_t i = 0; i < vec.size(); ++i) {
            result[i] = mod_(scalar * vec[i], q);
        }
        return result;
    };

    i128 q = pms.q;
    i128 p = pms.p;
    i128 g = pms.g;
    matrix_double A = pms.A;
    matrix_double B = pms.B;
    matrix_double C = pms.C;
    matrix_i128 F = pms.F;
    matrix_i128 G_bar = pms.G_bar;
    matrix_i128 H_bar = pms.H_bar;
    matrix_i128 R_bar = pms.R_bar;
    vector_double x_plant = pms.x_plant_init;
    vector_i128 x_cont = pms.x_cont_init_scaled;

    // i128 rr = pms.r;
    // i128 ss = pms.s;
    // i128 L = pms.L;
    // i128 iter_ = pms.iter_;

    i128 from = 1;
    i128 to_inclusive = 100;
    auto knowledge_exps = pms.sample_knowledge_exponents(from, to_inclusive);
    i128 alpha_0 = knowledge_exps.at(0);
    i128 alpha_1 = knowledge_exps.at(1);
    i128 gamma_0 = knowledge_exps.at(2);
    i128 gamma_1 = knowledge_exps.at(3);
    i128 rho_0 = knowledge_exps.at(4);
    i128 rho_1 = knowledge_exps.at(5);
    i128 m = H_bar.n_rows();
    i128 n = F.n_rows();
    auto verification_vectors = pms.sample_verification_vectors(m, n, from, to_inclusive);
    vector_i128 r_0 = verification_vectors.at(0);
    vector_i128 r_1 = verification_vectors.at(1);
    vector_i128 s = verification_vectors.at(2);

    const i128 N = 4;
    vector_i128 sk = sample_secret_key(N);
    const i128 d = 4;
    double log2q = std::log2(static_cast<double>(FIELD_MODULUS));
    int power = static_cast<int>(std::ceil(log2q / static_cast<double>(d)));
    i128 v = static_cast<i128>(1) << power;

    Encryptor enc(v, d, N, q, sk);
    i128 n_rlwes = 2 * d;
    i128 n_polys = 2;
    i128 n_coeffs = N;
    rgsw_mat F_ctx(F.n_rows(), F.n_cols(), n_rlwes, n_polys, n_coeffs);
    F_ctx = enc.encrypt_rgsw_mat(F);
    rgsw_mat G_bar_ctx(G_bar.n_rows(), G_bar.n_cols(), n_rlwes, n_polys, n_coeffs);
    G_bar_ctx = enc.encrypt_rgsw_mat(G_bar);
    rgsw_mat R_bar_ctx(R_bar.n_rows(), R_bar.n_cols(), n_rlwes, n_polys, n_coeffs);
    R_bar_ctx = enc.encrypt_rgsw_mat(R_bar);
    rgsw_mat H_bar_ctx(H_bar.n_rows(), H_bar.n_cols(), n_rlwes, n_polys, n_coeffs);
    H_bar_ctx = enc.encrypt_rgsw_mat(H_bar);
    rlwe_vec x_cont_ctx(x_cont.size(), n_polys, n_coeffs);
    x_cont_ctx = enc.encrypt_rlwe_vec(x_cont);
    rlwe_vec x_cont_ctx_convolved(x_cont_ctx); // copy x_cont_ctx

    vector_i128 eval_pows = eval_poly_pows(2 * N, 42, q);
    std::tuple<eval_key, veri_key> keys = compute_eval_and_veri_keys(
        F_ctx, G_bar_ctx, R_bar_ctx, H_bar_ctx,
        r_0, r_1, s, rho_0, rho_1, alpha_0, alpha_1, gamma_0, gamma_1,
        // pms.g, d, q, pms.p, v, N, eval_pows, enc
        g, d, q, p, N, eval_pows, enc
    );

    auto x_cont_ctx_hashed = x_cont_ctx_convolved.get_hash(eval_pows);
    // auto vec_dot_prod = [](vector_i128 scalar_vec, rlwe_vec rv) -> rlwe {
    //     assert(scalar_vec.size() == rv.size());
    //     rlwe dot_prod(rv.n_polys(), rv.n_coeffs());
    //     for (const auto& [x, y] : std::views::zip(scalar_vec, rv)) {
    //         for (size_t i = 0; i < y.size(); ++i) {
    //             dot_prod.set(i, dot_prod.get(i) + y.get(i) * x);
    //         }
    //     }
    //     return dot_prod;
    // };
    auto vec_dot_prod = [](vector_i128 scalar_vec, hashed_rlwe_vec rv) -> hashed_rlwe {
        assert(scalar_vec.size() == rv.size());
        hashed_rlwe dot_prod(rv.n_hashed_polys());
        for (const auto& [x, y] : std::views::zip(scalar_vec, rv)) {
            for (size_t i = 0; i < y.size(); ++i) {
                dot_prod.set(i, dot_prod.get(i) + y.get(i) * x);
            }
        }
        return dot_prod;
    };
    auto rx_0 = vec_dot_prod(r_1, x_cont_ctx_hashed);
    auto g1 = rx_0.pow();
    Proof old_proof {};
    old_proof.g_1 = g1;

    for (size_t k = 0; k < pms.iter_; k++) {
        // Helper function for matrix-vector multiplication (mod q)
        auto mat_vec_mult = [](const matrix_double& mat, const vector_double& vec, i128 q) -> vector_i128 {
            assert(mat[0].size() == vec.size());
            vector_i128 result(mat.size());
            for (size_t i = 0; i < mat.size(); ++i) {
                i128 sum = 0;
                for (size_t j = 0; j < mat[i].size(); ++j) {
                    sum = mod_(sum + static_cast<i128>(mat[i][j] * vec[j]), q);
                }
                result[i] = sum;
            }
            return result;
        };

        vector_i128 y_out = mat_vec_mult(C, x_plant, q);
        // Define round_vec to round each element of a vector<double> to i128
        auto round_vec = [](const vector_i128& vec) -> vector_i128 {
            vector_i128 res(vec.size());
            for (size_t i = 0; i < vec.size(); ++i) {
                res[i] = static_cast<i128>(std::llround(static_cast<double>(vec[i])));
            }
            return res;
        };
        y_out = scalar_vec_mult(pms.L, round_vec(scalar_vec_mult(pms.r, y_out, q)), q);
        rlwe_vec y_out_ctx = enc.encrypt_rlwe_vec(y_out);
        
    }
}

int main() {
    // test();
    omp_set_nested(1);
    // omp_set_num_threads(1);
    // std::cout << "omp num threads: " << omp_get_max_threads() << std::endl;
    // std::cout << " omp num threads: " << omp_get_num_threads() << std::endl;
    // test_full();
    run_control_loop();

    return 0;
}

// DONE update either i128 to u128, or replace % with mod()
// TODO test encryptor funcs (perhaps set moduli small and compare with Python)