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

// Anonymous function for cyclic exponentiation: (g^e1)^e2 mod p
i128 cyclic_exp(i128 base, i128 exp, i128 q, i128 p) {
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

std::tuple<std::tuple<eval_key, eval_key>, veri_key> compute_eval_and_veri_keys(
    rgsw_mat F_ctx, rgsw_mat G_bar_ctx, rgsw_mat R_bar_ctx, rgsw_mat H_bar_ctx,
    vector_i128 r_0, vector_i128 r_1, vector_i128 s,
    i128 rho_0, i128 rho_1, i128 alpha_0, i128 alpha_1, i128 gamma_0, i128 gamma_1,
    i128 g, i128 d, i128 q, i128 p, i128 N, vector_i128 eval_pows, Encryptor enc
    // i128 g, i128 d, i128 q, i128 p, i128 v, i128 N, vector_i128 eval_pows, Encryptor enc
) {
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
        rgsw_vec hashed(rgsw_v.size(), rgsw_v.n_rlwes(), rgsw_v.n_polys(), rgsw_v.n_coeffs());
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
    // rF_0 = hash_rgsw_vector(rF_0, eval_pows, q);
    rF_0 = hash_rgsw_vector(rF_0, eval_pows, q);
    rF_1 = hash_rgsw_vector(rF_1, eval_pows, q);

    // Make rgsw_vec for r_0 and r_1
    rgsw_vec r_0_rgsw(r_0.size(), 2 * d, 2, N);
    rgsw_vec r_1_rgsw(r_1.size(), 2 * d, 2, N);
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

    rgsw_vec rF_0_r_1(rF_0.size());
    rgsw_vec rF_1_r_0(rF_1.size());
    for (size_t i = 0; i < rF_0.size(); ++i)
        rF_0_r_1.set(i, rF_0.get(i) - r_1_rgsw.get(i));
    for (size_t i = 0; i < rF_1.size(); ++i)
        rF_1_r_0.set(i, rF_1.get(i) - r_0_rgsw.get(i));

    // = convert_vec_to_cyclic_enc(g, rF_0_r_1, d, q, p);
    hashed_rgsw_vec grFr_0 = rF_0_r_1.get_hash(eval_pows).pow();
    // convert_vec_to_cyclic_enc(g, rF_1_r_0, d, q, p);
    hashed_rgsw_vec grFr_1 = rF_1_r_0.get_hash(eval_pows).pow();

    // convert_vec_to_cyclic_enc(g, rF_0_r_1 * alpha_1, d, q, p);
    hashed_rgsw_vec grFr_alpha_0 = (rF_0_r_1.get_hash(eval_pows) * alpha_1).pow();
    // convert_vec_to_cyclic_enc(g, rF_1_r_0 * alpha_0, d, q, p);
    hashed_rgsw_vec grFr_alpha_1 = (rF_1_r_0.get_hash(eval_pows) * alpha_0).pow();

    hashed_rgsw_vec rG_0 = vec_mat_mult(r_0, G_bar_ctx).get_hash(eval_pows);
    hashed_rgsw_vec rG_1 = vec_mat_mult(r_1, G_bar_ctx).get_hash(eval_pows);
    hashed_rgsw_vec rR_0 = vec_mat_mult(r_0, R_bar_ctx).get_hash(eval_pows);
    hashed_rgsw_vec rR_1 = vec_mat_mult(r_1, R_bar_ctx).get_hash(eval_pows);

    // hashed_rgsw_vec sH = vec_mat_mult(s, H_bar_ctx).get_hash(eval_pows);
    rgsw_vec sH = vec_mat_mult(s, H_bar_ctx);

    rgsw_vec sH_r_1(sH.size());
    rgsw_vec sH_r_0(sH.size());
    for (size_t i = 0; i < sH.size(); ++i)
        sH_r_1.set(i, sH.get(i) - r_1_rgsw.get(i));
    for (size_t i = 0; i < sH.size(); ++i)
        sH_r_0.set(i, sH.get(i) - r_0_rgsw.get(i));

    // rgsw_vec gsHr_0 = convert_vec_to_cyclic_enc(g, mod_rgsw_vec(sH_r_1, q), d, q, p);
    hashed_rgsw_vec gsHr_0 = sH_r_1.get_hash(eval_pows).pow();
    // rgsw_vec gsHr_1 = convert_vec_to_cyclic_enc(g, mod_rgsw_vec(sH_r_0, q), d, q, p);
    hashed_rgsw_vec gsHr_1 = sH_r_0.get_hash(eval_pows).pow();

    // rgsw_vec gsHr_gamma_0 = convert_vec_to_cyclic_enc(g, mod_rgsw_vec(sH_r_1 * gamma_1, q), d, q, p);
    hashed_rgsw_vec gsHr_gamma_0 = (sH_r_1.get_hash(eval_pows) * gamma_1).pow();
    // rgsw_vec gsHr_gamma_1 = convert_vec_to_cyclic_enc(g, mod_rgsw_vec(sH_r_0 * gamma_0, q), d, q, p);
    hashed_rgsw_vec gsHr_gamma_1 = (sH_r_0.get_hash(eval_pows) * gamma_0).pow();



    eval_key ek0 {
        gr_0,
        grFr_0,
        gsHr_0,
        gr_rho_0,
        grFr_alpha_0,
        gsHr_gamma_0
    };

    eval_key ek1 {
        gr_1,
        grFr_1,
        gsHr_1,
        gr_rho_1,
        grFr_alpha_1,
        gsHr_gamma_1
    };

    auto ek = std::make_tuple(ek0, ek1);

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

// Helper function for scalar-vector multiplication (mod q)
template<typename vector_T>
auto scalar_vec_mult (double scalar, const vector_T& vec) {
    vector_double result(vec.size());
    for (size_t i = 0; i < vec.size(); ++i) {
        result[i] = scalar * vec[i];
    }
    return result;
};

// Computes proof values and returns a tuple of 2-tuples of i128
Proof compute_proof(
    const eval_key& ek,
    const rlwe_vec& x_nega_,
    const rlwe_vec& x_,
    const rlwe_vec& x,
    i128 v, i128 d,
    const vector_i128& eval_pows
) {
    const vector_i128& gr = ek.gr;
    const hashed_rgsw_vec& grFr = ek.grFr;
    const hashed_rgsw_vec& gsHr = ek.gsHr;
    const vector_i128& gr_rho = ek.gr_rho;
    const hashed_rgsw_vec& grFr_alpha = ek.grFr_alpha;
    const hashed_rgsw_vec& gsHr_gamma = ek.gsHr_gamma;

    // pow_() raises each element of the vector to the power of each element of the hashed_rlwe_decomp_vec
    // and returns a hashed_rlwe
    auto pow_ = [](const vector_i128& vec, const hashed_rlwe_vec& rv) -> hashed_rlwe {
        assert(vec.size() == rv.size());
        hashed_rlwe res(rv.n_hashed_polys());
        res.set_coeffs_to_one();
        for (size_t i = 0; i < vec.size(); ++i)
            res = res.group_mult(rv.get(i).pow(vec[i]));
        return res;
    };
    // cyclic_vec_dot_product returns a 2-tuple of i128
    auto grx_ = pow_(gr, x_.get_hash(eval_pows)); // G_1
    auto grFrx = grFr.pow(x.decompose(v, d).get_hash(eval_pows)); // G_2
    auto gsHrx = gsHr.pow(x.decompose(v, d).get_hash(eval_pows)); // G_3
    auto gr_rho_x_ = pow_(gr_rho, x_.get_hash(eval_pows)); // G_1_
    auto grFr_alpha_x = grFr_alpha.pow(x.decompose(v, d).get_hash(eval_pows)); // G_2_
    auto gsHr_gamma_x = gsHr_gamma.pow(x.decompose(v, d).get_hash(eval_pows)); // G_3_
    auto g_1 = pow_(gr, x_nega_.get_hash(eval_pows)); // g_1

    return Proof {
        grx_,
        grFrx,
        gsHrx,
        gr_rho_x_,
        grFr_alpha_x,
        gsHr_gamma_x,
        g_1,
    };
}

void verify_with_lin_and_dyn_checks(
    const veri_key& vk, const Proof& proof, const Proof& old_proof, i128 k,
    const rlwe_vec& y, const hashed_rlwe_vec& u, const rlwe_vec& u_reenc,
    i128 q, i128 p, const vector_i128& eval_pows
) {
    // Unpack veri_key
    const vector_i128& s = vk.s;
    const hashed_rgsw_vec& rG_0 = vk.rG_0;
    const hashed_rgsw_vec& rG_1 = vk.rG_1;
    const hashed_rgsw_vec& rR_0 = vk.rR_0;
    const hashed_rgsw_vec& rR_1 = vk.rR_1;
    i128 rho_0 = vk.rho_0;
    i128 rho_1 = vk.rho_1;
    i128 alpha_0 = vk.alpha_0;
    i128 alpha_1 = vk.alpha_1;
    i128 gamma_0 = vk.gamma_0;
    i128 gamma_1 = vk.gamma_1;

    // Unpack proofs
    const hashed_rlwe& g_1 = old_proof.g_1;
    const hashed_rlwe& G_1 = proof.grx_;
    const hashed_rlwe& G_2 = proof.grFrx;
    const hashed_rlwe& G_3 = proof.gsHrx;
    const hashed_rlwe& G_1_ = proof.gr_rho_x_;
    const hashed_rlwe& G_2_ = proof.grFr_alpha_x;
    const hashed_rlwe& G_3_ = proof.gsHr_gamma_x;

    // Select parameters based on k
    i128 rho, alpha, gamma;
    const hashed_rgsw_vec *rG, *rR;
    if (k % 2 == 0) {
        rho = rho_0;
        alpha = alpha_1;
        gamma = gamma_1;
        rG = &rG_0;
        rR = &rR_0;
    } else {
        rho = rho_1;
        alpha = alpha_0;
        gamma = gamma_0;
        rG = &rG_1;
        rR = &rR_1;
    }

    // Helper lambdas
    auto cyclic_exp_enc = [](const hashed_rlwe& hrlwe, i128 exp, i128 q, i128 p) -> hashed_rlwe {
        hashed_rlwe res(hrlwe.size());
        for (size_t i = 0; i < hrlwe.size(); ++i) {
            res.set(i, cyclic_exp(hrlwe.get(i), exp, q, p));
        }
        return res;
    };
    auto cyclic_mult_enc = [](const hashed_rlwe& a, const hashed_rlwe& b, i128 p) -> hashed_rlwe {
        assert(a.size() == b.size());
        hashed_rlwe res(a.size());
        for (size_t i = 0; i < a.size(); ++i) {
            res.set(i, mod_(a.get(i) * b.get(i), p));
        }
        return res;
    };
    auto vec_dot_prod = [](const vector_i128& vec, const hashed_rlwe_vec& hvec, i128 q) -> i128 {
        assert(vec.size() == hvec.size());
        i128 sum = 0;
        for (size_t i = 0; i < vec.size(); ++i) {
            for (size_t j = 0; j < hvec.get(i).size(); ++j) {
                sum = mod_(sum + mod_(vec[i] * hvec.get(i).get(j), q), q);
            }
        }
        return sum;
    };
    auto vec_dot_prod_enc = [](const hashed_rgsw_vec& rgsw_v, const rlwe_vec& rlwe_v, const vector_i128& eval_pows) -> hashed_rlwe {
        assert(rgsw_v.size() == rlwe_v.size());
        hashed_rlwe sum(N_POLYS_IN_RLWE);
        for (size_t i = 0; i < rgsw_v.size(); ++i) {
            // For each element, hash rlwe_v[i] and multiply by rgsw_v[i]
            auto hashed = rlwe_v.get(i).get_hash(eval_pows);
            for (size_t j = 0; j < hashed.size(); ++j) {
                sum = rgsw_v.get(i).get(j) * hashed.get(j);
            }
        }
        return sum;
    };

    // Linearity checks
    hashed_rlwe G_1_rho = cyclic_exp_enc(G_1, rho, q, p);
    hashed_rlwe G_2_alpha = cyclic_exp_enc(G_2, alpha, q, p);
    hashed_rlwe G_3_gamma = cyclic_exp_enc(G_3, gamma, q, p);
    for (size_t i = 0; i < 2; i++) {
        assert(G_1_rho.get(i) == G_1_.get(i));
        assert(G_2_alpha.get(i) == G_2_.get(i));
        assert(G_3_gamma.get(i) == G_3_.get(i));
    }

    // Dynamics check
    hashed_rlwe su = vec_dot_prod(s, u, q);
    hashed_rlwe gsu = su.pow();
    hashed_rlwe rhs_u = G_3.group_mult(g_1);
    for (size_t i = 0; i < 2; i++)
        assert(gsu.get(i) == rhs_u.get(i));

    hashed_rlwe rhs = G_2.group_mult(g_1);
    hashed_rlwe rGy = vec_dot_prod_enc(*rG, y, eval_pows);
    // hashed_rlwe grGy = cyclic_exp_enc(rhs, rGy, q, p);
    hashed_rlwe grGy = rhs.group_mult(rGy);
    rhs = cyclic_mult_enc(rhs, grGy, p);
    hashed_rlwe rRu = vec_dot_prod_enc(*rR, u_reenc, eval_pows);
    hashed_rlwe grRu = rRu.pow();
    rhs = rhs.group_mult(grRu);
    for (size_t i = 0; i < 2; i++)
        assert(G_1.get(i) == rhs.get(i));
}

void run_control_loop() {
    Params pms;
    using matrix_i128 = array2d<i128>;

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

    i128 rr = pms.r;
    i128 ss = pms.s;
    i128 L = pms.L;
    i128 iter_ = pms.iter_;

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
    rgsw_mat F_ctx = enc.encrypt_rgsw_mat(F);
    rgsw_mat G_bar_ctx = enc.encrypt_rgsw_mat(G_bar);
    rgsw_mat R_bar_ctx = enc.encrypt_rgsw_mat(R_bar);
    rgsw_mat H_bar_ctx = enc.encrypt_rgsw_mat(H_bar);
    rlwe_vec x_cont_ctx = enc.encrypt_rlwe_vec(x_cont);
    rlwe_vec x_cont_ctx_convolved(x_cont_ctx);

    i128 eval_point = 42;
    vector_i128 eval_pows = eval_poly_pows(2 * N, eval_point, q);
    auto keys = compute_eval_and_veri_keys(
        F_ctx, G_bar_ctx, R_bar_ctx, H_bar_ctx,
        r_0, r_1, s, rho_0, rho_1, alpha_0, alpha_1, gamma_0, gamma_1,
        g, d, q, p, N, eval_pows, enc
    );
    auto ek = std::get<0>(keys);
    veri_key vk = std::get<1>(keys);

    auto x_cont_ctx_hashed = x_cont_ctx_convolved.get_hash(eval_pows);

    auto vec_dot_prod = [](const vector_i128& scalar_vec, const hashed_rlwe_vec& rv) -> hashed_rlwe {
        assert(scalar_vec.size() == rv.size());
        hashed_rlwe dot_prod(rv.n_hashed_polys());
        for (size_t i = 0; i < scalar_vec.size(); ++i) {
            for (size_t j = 0; j < rv.get(i).size(); ++j) {
                dot_prod.set(j, mod_(dot_prod.get(j) + rv.get(i).get(j) * scalar_vec[i], FIELD_MODULUS));
            }
        }
        return dot_prod;
    };
    auto rx_0 = vec_dot_prod(r_1, x_cont_ctx_hashed);
    auto g1 = rx_0.pow();
    Proof old_proof {};
    old_proof.g_1 = g1;

    auto mat_vec_mult = [](const matrix_double& mat, const vector_double& vec) -> vector_double {
        assert(mat[0].size() == vec.size());
        vector_double result(mat.size());
        for (size_t i = 0; i < mat.size(); ++i) {
            double sum = 0;
            for (size_t j = 0; j < mat[i].size(); ++j) {
                sum += mat[i][j] * vec[j];
            }
            result[i] = sum;
        }
        return result;
    };

    auto round_and_mod_vec = [](const vector_double& vec) -> vector_i128 {
        vector_i128 res(vec.size());
        for (size_t i = 0; i < vec.size(); ++i) {
            i128 rounded = static_cast<i128>(std::round(vec[i]));
            i128 modded = mod_(rounded, FIELD_MODULUS);
            res[i] = modded;
        }
        return res;
    };

    for (size_t k = 0; k < iter_; k++) {
        // Plant: Compute output
        vector_double y_out = mat_vec_mult(C, x_plant);
        vector_i128 y_out_rounded_modded = round_and_mod_vec(scalar_vec_mult(L * rr, y_out));
        rlwe_vec y_out_ctx = enc.encrypt_rlwe_vec(y_out_rounded_modded); // P -> C

        // Controller: Compute output
        rlwe_vec u_out_ctx_convolved = H_bar_ctx.convolve(x_cont_ctx.decompose(v, d));
        hashed_rlwe_vec u_out_ctx_hashed = u_out_ctx_convolved.get_hash(eval_pows);
        rlwe_vec u_out_ctx = u_out_ctx_convolved.conv_to_nega(N); // C -> P

        // Plant: Update state (and re-encrypt u)
        vector_i128 u_out_ptx = enc.decrypt_rlwe_vec(u_out_ctx);
        vector_double u_in = scalar_vec_mult(1.0 / (rr * ss * ss * L), u_out_ptx);
        vector_i128 u_in_rounded_modded = round_and_mod_vec(scalar_vec_mult(rr * L, u_in));
        rlwe_vec u_reenc_ctx = enc.encrypt_rlwe_vec(u_in_rounded_modded); // XXX P->C: Plant re-encrypts u
        x_plant = mat_vec_mult(A, x_plant);

        // Controller: Update state
        rlwe_vec x_cont_old_ctx = x_cont_ctx;
        rlwe_vec x_cont_ctx_convolved =
            F_ctx.convolve(x_cont_ctx.decompose(v, d))
            + G_bar_ctx.convolve(y_out_ctx.decompose(v, d))
            + R_bar_ctx.convolve(u_reenc_ctx.decompose(v, d));
        x_cont_ctx = x_cont_ctx_convolved.conv_to_nega(N);

        // Controller: Prove
        eval_key& ek_i = std::get<0>(ek); // TODO check allowable &const variations
        if (k % 2 == 1)
            ek_i = std::get<1>(ek);
        Proof proof = compute_proof(ek_i, x_cont_ctx, x_cont_ctx_convolved, x_cont_old_ctx, v, d, eval_pows);

        // Plant: Verify
        verify_with_lin_and_dyn_checks(vk, proof, old_proof, k, y_out_ctx, u_out_ctx_hashed, u_reenc_ctx, q, p, eval_pows);
        old_proof = proof;
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