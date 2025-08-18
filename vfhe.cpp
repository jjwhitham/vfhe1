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
#include <iomanip>


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
    const rgsw_mat& F_ctx, const rgsw_mat& G_bar_ctx, const rgsw_mat& R_bar_ctx, const rgsw_mat& H_bar_ctx,
    const vector_i128& r_0, const vector_i128& r_1, const vector_i128& s,
    i128 rho_0, i128 rho_1, i128 alpha_0, i128 alpha_1, i128 gamma_0, i128 gamma_1,
    i128 d, i128 q, i128 N, const vector_i128& eval_pows, const Encryptor& enc
) {
    // Converts a vector of i128s v to a vector of g^v[i] mod p
    auto convert_vec_to_cyclic_a = [&](const vector_i128& v, const vector_i128& eval_pows) -> std::vector<hashed_a_poly> {
        std::vector<hashed_a_poly> res(v.size());
        // HACK ? halving eval_pows size
        i128 N = eval_pows.size() / 2;
        N = 2 * N - 1;
        hashed_a_poly eval_pows_a(N);
        for (size_t i = 0; i < N; i++)
            eval_pows_a.set(i, eval_pows.at(i));
        for (size_t i = 0; i < v.size(); ++i) {
            // mult v[i] by eval_pows
            hashed_a_poly hash_a = eval_pows_a * v[i];
            res.at(i) = hash_a.pow();
        }
        return res;
    };

    std::vector<hashed_a_poly> gr_0 = convert_vec_to_cyclic_a(r_0, eval_pows);
    std::vector<hashed_a_poly> gr_1 = convert_vec_to_cyclic_a(r_1, eval_pows);
    std::vector<hashed_a_poly> gr_rho_0 = convert_vec_to_cyclic_a(scalar_vec_mult(rho_0, r_0, q), eval_pows); // Example for first element
    std::vector<hashed_a_poly> gr_rho_1 = convert_vec_to_cyclic_a(scalar_vec_mult(rho_1, r_1, q), eval_pows);

    // Anonymous function for vector-matrix multiplication: vector_i128 * rgsw_mat
    auto vec_mat_mult = [](const vector_i128& vec, const rgsw_mat& mat) -> rgsw_vec {
        size_t rows = mat.n_rows();
        size_t cols = mat.n_cols();
        assert(vec.size() == rows);
        rgsw_vec res(cols, mat.n_rlwes(), mat.n_polys(), mat.n_coeffs());
        for (size_t j = 0; j < cols; ++j) {
            rgsw sum(mat.n_rlwes(), mat.n_polys(), mat.n_coeffs());
            for (size_t i = 0; i < rows; ++i) {
                sum = sum + (mat.get(i, j) * vec[i]);
            }
            res.set(j, sum);
        }
        return res;
    };

    hashed_a_rgsw_vec rF_0 = vec_mat_mult(r_0, F_ctx).get_hash_a(eval_pows);
    hashed_a_rgsw_vec rF_1 = vec_mat_mult(r_1, F_ctx).get_hash_a(eval_pows);

    // Make rgsw_vec for r_0 and r_1
    i128 n_hashed_a_coeffs = eval_pows.size() / 2;
    hashed_a_rgsw_vec r_0_rgsw(r_0.size(), N_POLYS_IN_RLWE * d, N_POLYS_IN_RLWE, n_hashed_a_coeffs);
    hashed_a_rgsw_vec r_1_rgsw(r_1.size(), N_POLYS_IN_RLWE * d, N_POLYS_IN_RLWE, n_hashed_a_coeffs);
    for (size_t i = 0; i < r_0.size(); ++i) {
        poly p(N);
        p.set(0, r_0[i]);
        r_0_rgsw.set(i, enc.encode_rgsw(p).get_hash_a(eval_pows));
    }
    for (size_t i = 0; i < r_1.size(); ++i) {
        poly p(N);
        p.set(0, r_1[i]);
        r_1_rgsw.set(i, enc.encode_rgsw(p).get_hash_a(eval_pows));
    }
    assert(r_0_rgsw.size() == rF_0.size());
    // assert(r_0_rgsw[0].rows() == rF_0[0].rows()); // FIXME

    hashed_a_rgsw_vec rF_0_r_1(rF_0.size());
    hashed_a_rgsw_vec rF_1_r_0(rF_1.size());
    for (size_t i = 0; i < rF_0.size(); ++i)
        rF_0_r_1.set(i, rF_0.get(i) - r_1_rgsw.get(i));
    for (size_t i = 0; i < rF_1.size(); ++i)
        rF_1_r_0.set(i, rF_1.get(i) - r_0_rgsw.get(i));

    hashed_a_rgsw_vec grFr_0 = rF_0_r_1.pow();
    hashed_a_rgsw_vec grFr_1 = rF_1_r_0.pow();
    hashed_a_rgsw_vec grFr_alpha_0 = (rF_0_r_1 * alpha_1).pow();
    hashed_a_rgsw_vec grFr_alpha_1 = (rF_1_r_0 * alpha_0).pow();

    hashed_rgsw_vec rG_0 = vec_mat_mult(r_0, G_bar_ctx).get_hash(eval_pows);
    hashed_rgsw_vec rG_1 = vec_mat_mult(r_1, G_bar_ctx).get_hash(eval_pows);
    hashed_rgsw_vec rR_0 = vec_mat_mult(r_0, R_bar_ctx).get_hash(eval_pows);
    hashed_rgsw_vec rR_1 = vec_mat_mult(r_1, R_bar_ctx).get_hash(eval_pows);

    hashed_a_rgsw_vec sH = vec_mat_mult(s, H_bar_ctx).get_hash_a(eval_pows);

    hashed_a_rgsw_vec sH_r_1(sH.size());
    hashed_a_rgsw_vec sH_r_0(sH.size());
    for (size_t i = 0; i < sH.size(); ++i)
        sH_r_1.set(i, sH.get(i) - r_1_rgsw.get(i));
    for (size_t i = 0; i < sH.size(); ++i)
        sH_r_0.set(i, sH.get(i) - r_0_rgsw.get(i));

    hashed_a_rgsw_vec gsHr_0 = sH_r_1.pow();
    hashed_a_rgsw_vec gsHr_1 = sH_r_0.pow();
    hashed_a_rgsw_vec gsHr_gamma_0 = (sH_r_1 * gamma_1).pow();
    hashed_a_rgsw_vec gsHr_gamma_1 = (sH_r_0 * gamma_0).pow();

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
    i128 v, i128 d
) {
    const std::vector<hashed_a_poly>& gr = ek.gr;
    const hashed_a_rgsw_vec& grFr = ek.grFr;
    const hashed_a_rgsw_vec& gsHr = ek.gsHr;
    const std::vector<hashed_a_poly>& gr_rho = ek.gr_rho;
    const hashed_a_rgsw_vec& grFr_alpha = ek.grFr_alpha;
    const hashed_a_rgsw_vec& gsHr_gamma = ek.gsHr_gamma;

    // pow_() raises each element of the vector to the power of each element of the hashed_rlwe_decomp_vec
    // and returns a hashed_rlwe
    auto pow_ = [](const std::vector<hashed_a_poly>& vec, const rlwe_vec& rv) -> hashed_rlwe {
        assert(vec.size() == rv.size());
        hashed_rlwe res(rv.n_polys());
        res.set_coeffs_to_one();
        for (size_t i = 0; i < vec.size(); i++) {
            hashed_a_poly hash_a_poly = vec.at(i);
            rlwe rl = rv.get(i);
            // TODO make refs?
            poly p1 = rl.get(0);
            poly p2 = rl.get(1);
            hashed_rlwe hash_rl(N_POLYS_IN_RLWE);
            hash_rl.set(0, hash_a_poly.get_hash_sec(p1));
            hash_rl.set(1, hash_a_poly.get_hash_sec(p2));
            res = res.group_mult(hash_rl);
        }
        return res;
    };

    auto grx_ = pow_(gr, x_); // G_1
    auto grFrx = grFr.get_hash_sec(x.decompose(v, d)); // G_2
    auto gsHrx = gsHr.get_hash_sec(x.decompose(v, d)); // G_3
    auto gr_rho_x_ = pow_(gr_rho, x_); // G_1_
    auto grFr_alpha_x = grFr_alpha.get_hash_sec(x.decompose(v, d)); // G_2_
    auto gsHr_gamma_x = gsHr_gamma.get_hash_sec(x.decompose(v, d)); // G_3_
    auto g_1 = pow_(gr, x_nega_); // g_1

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
    i128 v, i128 d, const vector_i128& eval_pows
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

    auto vec_dot_prod = [](const vector_i128& vec, const hashed_rlwe_vec& hvec) -> hashed_rlwe {
        assert(vec.size() == hvec.size());
        hashed_rlwe sum(N_POLYS_IN_RLWE);
        for (size_t i = 0; i < vec.size(); ++i) {
            sum = sum + hvec.get(i) * vec[i];
        }
        return sum;
    };
    auto vec_dot_prod_enc = [](const hashed_rgsw_vec& rgsw_v, const rlwe_vec& rlwe_v, const vector_i128& eval_pows, i128 v, i128 d) -> hashed_rlwe {
        assert(rgsw_v.size() == rlwe_v.size());
        hashed_rlwe sum(N_POLYS_IN_RLWE);
        for (size_t i = 0; i < rgsw_v.size(); ++i) {
            // For each element, hash rlwe_v[i] and multiply by rgsw_v[i]
            auto hashed = rlwe_v.get(i).decompose(v, d).get_hash(eval_pows);
            for (size_t j = 0; j < hashed.size(); ++j) {
                sum = sum + rgsw_v.get(i).get(j) * hashed.get(j);
            }
        }
        return sum;
    };

    // Linearity checks
    // hashed_rlwe G_1_rho = cyclic_exp_enc(G_1, rho, q, p);
    // TODO better interface self^other, with other as i128,
    // or better interface for i128^i128 (static member func?)
    DEBUG(std::cout << "G_1: ";)
    DEBUG(G_1.print();)
    DEBUG(std::cout << "\n";)
    DEBUG(std::cout << "G_1_: ";)
    DEBUG(G_1_.print();)
    DEBUG(std::cout << "\n";)
    for (size_t i = 0; i < 2; i++) {
        assert(G_1.pow_(G_1.get(i), rho) == G_1_.get(i));
        assert(G_2.pow_(G_2.get(i), alpha) == G_2_.get(i));
        assert(G_3.pow_(G_3.get(i), gamma) == G_3_.get(i));
    }

    // Dynamics check
    hashed_rlwe su = vec_dot_prod(s, u);
    hashed_rlwe gsu = su.pow();
    hashed_rlwe rhs_u = G_3.group_mult(g_1);
    for (size_t i = 0; i < 2; i++)
        assert(gsu.get(i) == rhs_u.get(i));

    hashed_rlwe rhs = G_2.group_mult(g_1);
    hashed_rlwe rGy = vec_dot_prod_enc(*rG, y, eval_pows, v, d);
    hashed_rlwe grGy = rGy.pow();
    rhs = rhs.group_mult(grGy);
    hashed_rlwe rRu = vec_dot_prod_enc(*rR, u_reenc, eval_pows, v, d);
    hashed_rlwe grRu = rRu.pow();
    rhs = rhs.group_mult(grRu);
    for (size_t i = 0; i < 2; i++)
        assert(G_1.get(i) == rhs.get(i));
}

void run_control_loop() {
    Params pms;
    DEBUG(pms.print();)
    DEBUG(std::cout << "*** START run_control_loop ***\n";)
    using matrix_i128 = array2d<i128>;

    i128 q = pms.q;
    matrix_double A = pms.A;
    matrix_double B = pms.B;
    matrix_double C = pms.C;
    matrix_i128 F = pms.F;
    matrix_i128 G_bar = pms.G_bar;
    matrix_i128 H_bar = pms.H_bar;
    matrix_i128 R_bar = pms.R_bar;
    vector_double x_plant = pms.x_plant_init;
    vector_i128 x_cont = pms.x_cont_init_scaled;
    DEBUG(std::cout << "C:\n";)
    DEBUG(for (auto& x : C) print_vector_double(x);)
    DEBUG(std::cout << "x_plant:\n";)
    DEBUG(print_vector_double(x_plant);)
    i128 rr = pms.r;
    i128 ss = pms.s;
    i128 L = pms.L;
    i128 iter_ = pms.iter_;
    DEBUG(iter_ = 1000;)

    i128 from = 1;
    i128 to_inclusive = q - 1;
    // TODO Params should sample
    auto knowledge_exps = pms.sample_knowledge_exponents(from, to_inclusive);
    knowledge_exps = vector_i128({1, 1, 3, 3, 2, 3});
    // DEBUG(std::cout << "knowledge_exps: ";)
    // DEBUG(print_vector_i128(knowledge_exps);)
    i128 alpha_0 = knowledge_exps.at(0);
    i128 alpha_1 = knowledge_exps.at(1);
    i128 gamma_0 = knowledge_exps.at(2);
    i128 gamma_1 = knowledge_exps.at(3);
    i128 rho_0 = knowledge_exps.at(4);
    i128 rho_1 = knowledge_exps.at(5);
    i128 m = H_bar.n_rows();
    i128 n = F.n_rows();
    // TODO Params should sample
    // r_0: 1, 2, 4, 3, 11
    // r_1: 2, 4, 1, 4, 2,
    // s: 1, 3,
    auto verification_vectors = pms.sample_verification_vectors(m, n, from, to_inclusive);
    vector_i128 r_0 = verification_vectors.at(0);
    r_0 = vector_i128({1, 2, 4, 3, 1});
    // DEBUG(std::cout << "r_0: ";)
    // DEBUG(print_vector_i128(r_0);)
    vector_i128 r_1 = verification_vectors.at(1);
    r_1 = vector_i128({2, 4, 1, 4, 2});
    // DEBUG(std::cout << "r_1: ";)
    // DEBUG(print_vector_i128(r_1);)
    vector_i128 s = verification_vectors.at(2);
    s = vector_i128({1, 3});
    // DEBUG(std::cout << "s: ";)
    // DEBUG(print_vector_i128(s);)

    const i128 N = N_;
    DEBUG1(const i128 N = 1 << 2;)
    // DEBUG(std::cout << "N_THREADS: ";)
    // DEBUG(std::cout << N_THREADS << "\n" ;)
    // DEBUG(std::cout << "N: ";)
    // DEBUG(std::cout << print_to_string_i128(N) << "\n" ;)
    // TODO move to Encryptor
    vector_i128 sk = sample_secret_key(N);
    DEBUG1(sk = vector_i128({0, 0});)
    // DEBUG(std::cout << "sk: ";)
    // DEBUG(print_vector_i128(sk);)
    // TODO move to Params
    const i128 d = 4;
    // const i128 d = 2;
    // DEBUG(std::cout << "d: ";)
    // DEBUG(std::cout << print_to_string_i128(d) << "\n" ;)
    // TODO move to Params
    double log2q = std::log2(static_cast<double>(FIELD_MODULUS));
    // DEBUG(std::cout << "log2q: ";)
    // DEBUG(std::cout << log2q << "\n" ;)
    // TODO move to Params
    int power = static_cast<int>(std::ceil(log2q / static_cast<double>(d)));
    i128 v = static_cast<i128>(1) << power;
    // DEBUG(std::cout << "v: ";)
    // DEBUG(std::cout << print_to_string_i128(v) << "\n" ;)

    Encryptor enc(v, d, N, q, sk);
    rgsw_mat F_ctx = enc.encrypt_rgsw_mat(F);
    // DEBUG(std::cout << "F_ctx.get(0, 0):\n";)
    // DEBUG(F_ctx.get(0, 0).print();)
    rgsw_mat G_bar_ctx = enc.encrypt_rgsw_mat(G_bar);
    // DEBUG(std::cout << "G_bar_ctx[0][0]:\n";)
    // DEBUG(G_bar_ctx.get(0, 0).print();)
    rgsw_mat R_bar_ctx = enc.encrypt_rgsw_mat(R_bar);
    rgsw_mat H_bar_ctx = enc.encrypt_rgsw_mat(H_bar);
    // DEBUG(std::cout << "H_bar_ctx.get(0, 0):\n";)
    // DEBUG(H_bar_ctx.get(0, 0).print();)
    rlwe_vec x_cont_ctx = enc.encrypt_rlwe_vec(x_cont);
    DEBUG(std::cout << "x_cont_ctx:\n";)
    DEBUG(x_cont_ctx.print();)
    rlwe_vec x_cont_ctx_convolved(x_cont_ctx);
    DEBUG(std::cout << "x_cont_ctx_convolved:\n";)
    DEBUG(x_cont_ctx_convolved.print();)


    i128 eval_point = 42;
    vector_i128 eval_pows = eval_poly_pows(2 * N, eval_point, q);
    auto keys = compute_eval_and_veri_keys(
        F_ctx, G_bar_ctx, R_bar_ctx, H_bar_ctx,
        r_0, r_1, s, rho_0, rho_1, alpha_0, alpha_1, gamma_0, gamma_1,
        d, q, N, eval_pows, enc
    );
    auto ek = std::get<0>(keys);
    veri_key vk = std::get<1>(keys);

    auto x_cont_ctx_hashed = x_cont_ctx.get_hash(eval_pows);

    auto vec_dot_prod = [](const vector_i128& vec, const hashed_rlwe_vec& hvec) -> hashed_rlwe {
        assert(vec.size() == hvec.size());
        hashed_rlwe sum(N_POLYS_IN_RLWE);
        for (size_t i = 0; i < vec.size(); ++i) {
            sum = sum + hvec.get(i) * vec[i];
        }
        return sum;
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

    TIMING(times_counts.iter_ = iter_;)
    for (size_t k = 0; k < iter_; k++) {
        DEBUG(std::cout << "\n*** START k: " << k << ", run_control_loop ***\n";)
        // Plant: Compute output
        TIMING(auto start_plant = std::chrono::high_resolution_clock::now();)

        vector_double y_out = mat_vec_mult(C, x_plant);
        DEBUG(std::cout << "y_out:\n";)
        DEBUG(print_vector_double(y_out);)
        vector_double y_out_scaled = scalar_vec_mult(rr, y_out);
        DEBUG(std::cout << "y_out_scaled:\n";)
        DEBUG(print_vector_double(y_out_scaled);)
        vector_i128 y_out_rounded_modded = round_and_mod_vec(y_out_scaled);
        DEBUG(std::cout << "y_out_scaled_modded:\n";)
        DEBUG(print_vector_i128(y_out_rounded_modded);)
        y_out_scaled = scalar_vec_mult(L, y_out_rounded_modded);
        DEBUG(std::cout << "y_out_scaled:\n";)
        DEBUG(print_vector_double(y_out_scaled);)
        y_out_rounded_modded = round_and_mod_vec(y_out_scaled);
        DEBUG(std::cout << "y_out_scaled_modded:\n";)
        DEBUG(print_vector_i128(y_out_rounded_modded);)
        rlwe_vec y_out_ctx = enc.encrypt_rlwe_vec(y_out_rounded_modded); // P -> C

        #ifdef TIMING_ON
            auto end_plant = std::chrono::high_resolution_clock::now();
            times_counts.elapsed_plant += end_plant - start_plant;
        #endif

        // Controller: Compute output
        rlwe_decomp_vec x_cont_ctx_decomp = x_cont_ctx.decompose(v, d);
        // DEBUG(std::cout << "x_cont_ctx_decomp:\n";)
        // DEBUG(x_cont_ctx_decomp.print();)
        rlwe_vec u_out_ctx_convolved = H_bar_ctx.convolve(x_cont_ctx_decomp); // C -> P
        // DEBUG(std::cout << "u_out_ctx_convolved:\n";)
        // DEBUG(u_out_ctx_convolved.print();)

        // Plant: Update state (and re-encrypt u)
        TIMING(start_plant = std::chrono::high_resolution_clock::now();)

        hashed_rlwe_vec u_out_ctx_hashed = u_out_ctx_convolved.get_hash(eval_pows);
        // DEBUG(std::cout << "u_out_ctx_hashed:\n";)
        // DEBUG(u_out_ctx_hashed.print();)
        rlwe_vec u_out_ctx = u_out_ctx_convolved.conv_to_nega(N);
        // DEBUG(std::cout << "u_out_ctx:\n";)
        // DEBUG(u_out_ctx.print();)
        vector_i128 u_out_ptx = enc.decrypt_rlwe_vec(u_out_ctx);
        // DEBUG(std::cout << "u_out_ptx:\n";)
        // DEBUG(print_vector_i128(u_out_ptx);)
        vector_double u_in = scalar_vec_mult(1.0 / (rr * ss * ss * L), u_out_ptx);
        DEBUG(std::cout << "u_in:\n";)
        DEBUG(print_vector_double(u_in);)
        vector_i128 u_in_rounded_modded = round_and_mod_vec(scalar_vec_mult(rr, u_in));
        u_in_rounded_modded = round_and_mod_vec(scalar_vec_mult(L, u_in));
        rlwe_vec u_reenc_ctx = enc.encrypt_rlwe_vec(u_in_rounded_modded); // XXX P->C: Plant re-encrypts u
        x_plant = mat_vec_mult(A, x_plant);
        vector_double B_u = mat_vec_mult(B, u_in);
        for (size_t i = 0; i < B_u.size(); i++)
            x_plant.at(i) += B_u.at(i);
        DEBUG(std::cout << "x_plant:\n";)
        DEBUG(print_vector_double(x_plant);)

        #ifdef TIMING_ON
            end_plant = std::chrono::high_resolution_clock::now();
            times_counts.elapsed_plant += end_plant - start_plant;
        #endif

        // Controller: Update state
        TIMING(auto start_controller = std::chrono::high_resolution_clock::now();)

        rlwe_vec x_cont_old_ctx = x_cont_ctx;
        DEBUG(std::cout << "x_cont_old_ctx:\n";)
        DEBUG(x_cont_old_ctx.print();)
        rlwe_vec F_x = F_ctx.convolve(x_cont_ctx.decompose(v, d));
        // DEBUG(std::cout << "F_x:\n";)
        // DEBUG(F_x.print();)
        rlwe_decomp_vec y_out_ctx_decomped = y_out_ctx.decompose(v, d);
        rlwe_vec G_y = G_bar_ctx.convolve(y_out_ctx_decomped);
        // DEBUG(std::cout << "y_out_ctx_decomped:\n";)
        // DEBUG(y_out_ctx_decomped.print();)
        // DEBUG(std::cout << "G_bar_ctx:";)
        // DEBUG(G_bar_ctx.print();)
        // DEBUG(std::cout << "y_out_ctx:\n";)
        // DEBUG(y_out_ctx.print();)
        // DEBUG(std::cout << "G_y:\n";)
        // DEBUG(G_y.print();)
        rlwe_vec R_u = R_bar_ctx.convolve(u_reenc_ctx.decompose(v, d));
        // DEBUG(std::cout << "u_reenc_ctx:\n";)
        // DEBUG(u_reenc_ctx.print();)
        // DEBUG(std::cout << "R_u:\n";)
        // DEBUG(R_u.print();)
        rlwe_vec x_cont_ctx_convolved = F_x + G_y + R_u;
        x_cont_ctx = x_cont_ctx_convolved.conv_to_nega(N);
        DEBUG(std::cout << "x_cont_ctx_convolved:\n";)
        DEBUG(x_cont_ctx_convolved.print();)
        DEBUG(std::cout << "x_cont_ctx:\n";)
        DEBUG(x_cont_ctx.print();)

        #ifdef TIMING_ON
            auto end_controller = std::chrono::high_resolution_clock::now();
            times_counts.elapsed_controller += end_controller - start_controller;
        #endif

        // Controller: Prove
        #ifdef TIMING_ON
            start_controller = std::chrono::high_resolution_clock::now();
            auto start_proof = std::chrono::high_resolution_clock::now();
        #endif

        const eval_key& ek_i = (k % 2 == 0) ? std::get<0>(ek) : std::get<1>(ek);
        Proof proof = compute_proof(ek_i, x_cont_ctx, x_cont_ctx_convolved, x_cont_old_ctx, v, d); // C -> P

#ifdef TIMING_ON
            auto end_proof = std::chrono::high_resolution_clock::now();
            times_counts.elapsed_proof += end_proof - start_proof;
            end_controller = std::chrono::high_resolution_clock::now();
            times_counts.elapsed_controller += end_controller - start_controller;
        #endif


        // Plant: Verify
        #ifdef TIMING_ON
            start_plant = std::chrono::high_resolution_clock::now();
            auto start_verify = std::chrono::high_resolution_clock::now();
        #endif

        verify_with_lin_and_dyn_checks(vk, proof, old_proof, k, y_out_ctx, u_out_ctx_hashed, u_reenc_ctx, v, d, eval_pows);

        #ifdef TIMING_ON
            auto end_verify = std::chrono::high_resolution_clock::now();
            times_counts.elapsed_verify += end_verify - start_verify;
            end_plant = std::chrono::high_resolution_clock::now();
            times_counts.elapsed_plant += end_plant - start_plant;
        #endif

        old_proof = proof;
        DEBUG(std::cout << "\n\n*** END k: " << k << ", run_control_loop ***\n\n";)
    }
}

void print_times_and_counts() {
    int iter_ = times_counts.iter_;
    int n_decimals = 0;
    std::cout << std::fixed << std::setprecision(n_decimals);
    std::cout << "Times:\n";
    std::cout << "Proof (per loop): ";
    std::cout << times_counts.elapsed_proof.count() / iter_ << "\n";
    std::cout << "Controller (per loop): ";
    std::cout << times_counts.elapsed_controller.count() / iter_ << "\n";
    std::cout << "Verify (per loop): ";
    std::cout << times_counts.elapsed_verify.count() / iter_ << "\n";
    std::cout << "Plant (per loop): ";
    std::cout << times_counts.elapsed_plant.count() / iter_ << "\n";
    std::cout << "Total Elapsed time: ";
    std::cout << times_counts.elapsed_total.count() << "\n";
    std::cout << "Total Elapsed time (per loop): ";
    std::cout << times_counts.elapsed_total.count() / iter_ << "\n";

    std::cout << "\n";
    double convolve = 0.0;
    for (size_t i = 0; i < N_THREADS; i++)
        convolve += times_counts.convolve[i].count();
    std::cout << "Convolve (per loop): ";
    std::cout  << convolve / iter_ << "\n";
    std::cout << "Convolve (per loop, per thread): ";
    std::cout << convolve / N_THREADS / iter_ << "\n";
    std::cout << "Number of calls:\n";
    int convolve_calls = 0;
    for (size_t i = 0; i < N_THREADS; i++)
        convolve_calls += times_counts.calls_convolve[i];
    std::cout << "  Convolve (per loop): ";
    std::cout << convolve_calls / iter_ << "\n";
}

int main() {
    // omp_set_nested(1);
    TIMING(auto start = std::chrono::high_resolution_clock::now();)

    run_control_loop();

    #ifdef TIMING_ON
        auto end = std::chrono::high_resolution_clock::now();
        times_counts.elapsed_total = end - start;
        print_times_and_counts();
    #endif
    return 0;
}

// DONE timings
//  TODO csv output for timings
// TODO 128bit atcoder
// TODO fix verify bug
// TODO multiprecision library
// TODO add secure hash
// TODO increase group size
// TODO implement packing?
// TODO keep polys in NTT form?
    // TODO at what points in the code?
// TODO what happens when trying to use a cyclic group
//      when Q=q1*...*qn? Refer to RSA which has composite modulus
