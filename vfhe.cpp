#include "shared.h"
#include "enc.h"
#include "vfhe.h"
#include "ntt.h"
#include <ranges>
#include <iomanip>
#include "/Users/jw/Projects/mcl/include/mcl/bn.hpp"

using namespace mcl::bn;

std::tuple<std::tuple<eval_key, eval_key>, veri_key, check_key> compute_eval_and_veri_keys(
    const rgsw_mat& F_ctx, const rgsw_mat& G_bar_ctx, const rgsw_mat& R_bar_ctx, const rgsw_mat& H_bar_ctx,
    const veri_vec& r_0, const veri_vec& r_1, const veri_vec& s,
    bigz rho_0, bigz rho_1, bigz alpha_0, bigz alpha_1, bigz gamma_0, bigz gamma_1,
    size_t N, const vector_bigz& eval_pows, const Encryptor& enc
) {
    // TODO can try hash rgsw's first, then subtract r_0/r_1
    ASSERT(r_0.size() == r_1.size());
    flat_rgsw_vec rF_0 = (r_0 * F_ctx);
    flat_rgsw_vec rF_1 = (r_1 * F_ctx);
    ASSERT(r_0.size() == rF_0.size());
    // ASSERT(rF_1.n_polys() == (size_t)(2 * d)); // XXX
    ASSERT(rF_1.n_flat_rgsws() == F_ctx.n_cols());
    flat_rgsw_vec r_0_rgsw(r_0.size());
    flat_rgsw_vec r_1_rgsw(r_1.size());
    poly p{N};
    p.set(0, bigz{1});
    for (size_t i = 0; i < r_0.size(); i++) {
        r_0_rgsw.set(i, enc.encode_flat_rgsw(p * r_0.get(i).get(0), p * r_0.get(i).get(1)));
        // XXX
        r_1_rgsw.set(i, enc.encode_flat_rgsw(p * r_1.get(i).get(0), p * r_1.get(i).get(1)));
    }
    // hashed_t_rgsw_vec rF_0_r_1 = (rF_0 - r_1_rgsw).get_hash_t(eval_pows);
    hashed_rgsw_vec rF_0_r_1 = (rF_0 - r_1_rgsw).get_hash(eval_pows);
    // hashed_t_rgsw_vec rF_1_r_0 = (rF_1 - r_0_rgsw).get_hash_t(eval_pows);
    hashed_rgsw_vec rF_1_r_0 = (rF_1 - r_0_rgsw).get_hash(eval_pows);
    hashed_t_rgsw_vec grFr_0 = rF_0_r_1.pow();
    hashed_t_rgsw_vec grFr_1 = rF_1_r_0.pow();
    hashed_t_rgsw_vec grFr_alpha_0 = (rF_0_r_1 * alpha_1).pow();
    hashed_t_rgsw_vec grFr_alpha_1 = (rF_1_r_0 * alpha_0).pow();

    flat_rgsw_vec sH = s * H_bar_ctx;
    // hashed_t_rgsw_vec sH_r_1 = (sH - r_1_rgsw).get_hash_t(eval_pows);
    hashed_rgsw_vec sH_r_1 = (sH - r_1_rgsw).get_hash(eval_pows);
    // hashed_t_rgsw_vec sH_r_0 = (sH - r_0_rgsw).get_hash_t(eval_pows);
    hashed_rgsw_vec sH_r_0 = (sH - r_0_rgsw).get_hash(eval_pows);
    hashed_t_rgsw_vec gsHr_0 = sH_r_1.pow();
    hashed_t_rgsw_vec gsHr_1 = sH_r_0.pow();
    hashed_t_rgsw_vec gsHr_gamma_0 = (sH_r_1 * gamma_1).pow();
    hashed_t_rgsw_vec gsHr_gamma_1 = (sH_r_0 * gamma_0).pow();


    using ht_veri_vec = hashed_t_veri_vec;
    ht_veri_vec gr_0 = r_0.pow(); // TODO change get_hash_t to copy
    ht_veri_vec gr_1 = r_1.pow();
    ht_veri_vec gr_rho_0 = (r_0 * rho_0).pow();
    ht_veri_vec gr_rho_1 = (r_1 * rho_1).pow();

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

    // comput veri key
    hashed_rgsw_vec rG_0 = (r_0 * G_bar_ctx).get_hash(eval_pows);
    hashed_rgsw_vec rG_1 = (r_1 * G_bar_ctx).get_hash(eval_pows);
    hashed_rgsw_vec rR_0 = (r_0 * R_bar_ctx).get_hash(eval_pows);
    hashed_rgsw_vec rR_1 = (r_1 * R_bar_ctx).get_hash(eval_pows);
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

    check_key ck {
        r_1_rgsw,
        r_0_rgsw
    };

    return std::make_tuple(ek, vk, ck);
}


// Helper function for scalar-vector multiplication (mod q)
template<typename vector_T>
auto scalar_vec_mult (double scalar, const vector_T& vec) {
    vector_double result(vec.size());
    for (size_t i = 0; i < vec.size(); i++) {
        result[i] = scalar * vec[i];
    }
    return result;
};

// Computes proof values and returns a tuple of 2-tuples of bigz
Proof compute_proof(
    const eval_key& ek,
    rlwe_vec& x_nega_lhs,
    rlwe_vec& x_hi,
    const rlwe_vec& x,
    const bigz& v, u32 d, u32 power,
    std::vector<G1>& eval_pows_g
) {
    const hashed_t_veri_vec& gr = ek.gr;
    const hashed_t_rgsw_vec& grFr = ek.grFr;
    const hashed_t_rgsw_vec& gsHr = ek.gsHr;
    const hashed_t_veri_vec& gr_rho = ek.gr_rho;
    const hashed_t_rgsw_vec& grFr_alpha = ek.grFr_alpha;
    const hashed_t_rgsw_vec& gsHr_gamma = ek.gsHr_gamma;

    rlwe_decomp_vec x_decomped = x.decompose(v, d, power);
    x_decomped.msm(eval_pows_g);
    x_hi.msm(eval_pows_g);
    x_nega_lhs.msm(eval_pows_g);

    auto grx_ = gr.get_hash_sec(x_nega_lhs); // XXX G_1
    auto grx_hi = gr.get_hash_sec(x_hi); // XXX G_1_hi
    auto grFrx = grFr.get_hash_sec(x_decomped); // G_2
    auto gsHrx = gsHr.get_hash_sec(x_decomped); // G_3
    auto gr_rho_x_ = gr_rho.get_hash_sec(x_nega_lhs); // XXX G_1_
    auto gr_rho_x_hi = gr_rho.get_hash_sec(x_hi); // XXX G_1_
    auto grFr_alpha_x = grFr_alpha.get_hash_sec(x_decomped); // G_2_
    auto gsHr_gamma_x = gsHr_gamma.get_hash_sec(x_decomped); // G_3_

    return Proof {
        grx_,
        grx_hi,
        grFrx,
        gsHrx,
        gr_rho_x_,
        gr_rho_x_hi,
        grFr_alpha_x,
        gsHr_gamma_x,
    };
}

void check_assert(bool outcome) {
    if (outcome == true)
        std::cout << "PASSED\n";
    else if (outcome == false)
        std::cout << "FAILED\n";
    else
        std::cout << "ELSE\n";
}

void verify_with_lin_and_dyn_checks(
    const veri_key& vk, const Proof& proof, const Proof& old_proof, size_t k,
    const rlwe_vec& y, const rlwe_vec& u_conv, const rlwe_vec& u_reenc,
    bigz v, u32 d, u32 power, const vector_bigz& eval_pows
) {
    // Unpack veri_key
    const veri_vec& s = vk.s;
    const hashed_rgsw_vec& rG_0 = vk.rG_0;
    const hashed_rgsw_vec& rG_1 = vk.rG_1;
    const hashed_rgsw_vec& rR_0 = vk.rR_0;
    const hashed_rgsw_vec& rR_1 = vk.rR_1;
    bigz rho_0 = vk.rho_0;
    bigz rho_1 = vk.rho_1;
    bigz alpha_0 = vk.alpha_0;
    bigz alpha_1 = vk.alpha_1;
    bigz gamma_0 = vk.gamma_0;
    bigz gamma_1 = vk.gamma_1;

    // Unpack proofs
    const GT& g_1 = old_proof.grx_;
    const GT& G_1 = proof.grx_;
    const GT& G_hi = proof.grx_hi;
    const GT& G_2 = proof.grFrx;
    const GT& G_3 = proof.gsHrx;
    const GT& G_1_ = proof.gr_rho_x_;
    const GT& G_hi_ = proof.gr_rho_x_hi;
    const GT& G_2_ = proof.grFr_alpha_x;
    const GT& G_3_ = proof.gsHr_gamma_x;

    // Select parameters based on k
    bigz rho, alpha, gamma;
    const hashed_rgsw_vec &rG = (k % 2 == 0) ? rG_0 : rG_1;
    const hashed_rgsw_vec &rR = (k % 2 == 0) ? rR_0 : rR_1;
    if (k % 2 == 0) {
        rho = rho_0;
        alpha = alpha_1;
        gamma = gamma_1;
    } else {
        rho = rho_1;
        alpha = alpha_0;
        gamma = gamma_0;
    }

    static bigz tN_1_ = eval_pows.at(N_) + 1;
    static GT lhs;
    static GT rhs;
    static GT rhs1;
    // Linearity checks
    // Dynamics checks: controller output
    // Dynamics checks: controller state update

    static std::chrono::duration<double, std::milli> times[4];
    // #pragma omp parallel sections
    {
        // #pragma omp section
        {
            TIMING(auto start = std::chrono::high_resolution_clock::now();)
            GT gsu = pow_t(genT, s.dot_prod(u_conv.get_hash(eval_pows)));
            GT rhs_u = G_3 * g_1;
            assert(gsu == rhs_u);
            TIMING(auto end = std::chrono::high_resolution_clock::now();)
            times[0] += end - start;
        }
        // #pragma omp section
        {
            TIMING(auto start = std::chrono::high_resolution_clock::now();)
            rhs = G_2 * pow_t(genT, rG.dot_prod(y.decompose(v, d, power).get_hash(eval_pows)));
            TIMING(auto end = std::chrono::high_resolution_clock::now();)
            times[1] += end - start;
        }
        // #pragma omp section
        {
            TIMING(auto start = std::chrono::high_resolution_clock::now();)
            rhs1 = g_1 * pow_t(genT, rR.dot_prod(u_reenc.decompose(v, d, power).get_hash(eval_pows)));
            TIMING(auto end = std::chrono::high_resolution_clock::now();)
            times[2] += end - start;
        }
        // #pragma omp section
        {
            TIMING(auto start = std::chrono::high_resolution_clock::now();)
            assert(pow_t(G_1, rho) == G_1_);
            assert(pow_t(G_hi, rho) == G_hi_);
            assert(pow_t(G_2, alpha) == G_2_);
            assert(pow_t(G_3, gamma) == G_3_);
            lhs = G_1 * pow_t(G_hi, tN_1_);
            TIMING(auto end = std::chrono::high_resolution_clock::now();)
            times[3] += end - start;
        }
    }
    assert(lhs == rhs * rhs1);
    for (size_t i = 0; i < 4; i++)
        std::cout << "veri_thread" << i << ": " << times[i].count() << "\n";
}


struct control_law_vars {
    std::vector<vector_double> y;
    std::vector<vector_double> u;
    std::vector<vector_double> x_plant;
    std::vector<vector_double> x_cont;
};

void print_times_and_counts(times_and_counts& timing) {
    int iter_ = timing.iter_;
    int n_decimals = 0;
    std::cout << std::fixed << std::setprecision(n_decimals);
    std::cout << "\nTimes:\n";
    std::cout << "Total: ";
    std::cout << timing.total.count() << "\n";
    std::cout << "Setup: ";
    std::cout << timing.total.count() - timing.loop.count() << "\n";
    std::cout << "Per Loop: ";
    std::cout << timing.loop.count() / iter_ << "\n";
    std::cout << "  Controller: ";
    std::cout << timing.controller.count() / iter_ << "\n";
    std::cout << "    Proof: ";
    std::cout << timing.proof.count() / iter_ << "\n";
    std::cout << "  Plant: ";
    std::cout << timing.plant.count() / iter_ << "\n";
    std::cout << "    Verify: ";
    std::cout << timing.verify.count() / iter_ << "\n";
    std::cout << "  MSMs setup: ";
    std::cout << timing.msm.count() / iter_ << "\n";
    std::cout << "  MSMs compute: ";
    std::cout << timing.msm1.count() / iter_ << "\n";
    std::cout << "  NTT: ";
    std::cout << timing.ntt.count() / iter_ << "\n";
    std::cout << "  iNTT: ";
    std::cout << timing.intt.count() / iter_ << "\n";
    std::cout << "  NTT1: ";
    std::cout << timing.ntt1.count() / iter_ << "\n";
    std::cout << "  iNTT1: ";
    std::cout << timing.intt1.count() / iter_ << "\n";
    std::cout << "\n";

    std::cout << "Number of calls (per loop):\n";

    std::cout << "  NTT: ";
    std::cout << timing.calls_ntt / iter_ << "\n";
    std::cout << "  iNTT: ";
    std::cout << timing.calls_intt / iter_ << "\n";
    std::cout << "  NTT1: ";
    std::cout << timing.calls_ntt1 / iter_ << "\n";
    std::cout << "  iNTT1: ";
    std::cout << timing.calls_intt1 / iter_ << "\n";
    std::cout << "  Conv-to-nega: ";
    std::cout << timing.calls_conv_to_nega / iter_ << "\n";
    std::cout << "  MSMs: ";
    std::cout << timing.calls_msm / iter_ << "\n\n";
}

vector_double mat_vec_mult(const matrix_double& mat, const vector_double& vec) {
    ASSERT(mat[0].size() == vec.size());
    vector_double result(mat.size());
    for (size_t i = 0; i < mat.size(); i++) {
        mpf_class sum = 0;
        for (size_t j = 0; j < mat[i].size(); j++) {
            sum += mat[i][j] * vec[j];
        }
        result[i] = sum;
    }
    return result;
}

void run_control_loop(control_law_vars& vars, times_and_counts& timing) {

    auto eval_poly_pows = [](size_t n, const bigz& eval_point_t, const bigz& q) -> vector_bigz {
        vector_bigz res(n);
        res.at(0) = 1;
        for (size_t i = 1; i < n - 1; i++) { // NOTE last elem should be 0 for conv
            res.at(i) = res.at(i - 1) * eval_point_t;
            res.at(i) = mod_(res.at(i), q);
        }
        assert(res.at(res.size() - 1) == 0);
        return res;
    };
    auto make_eval_pows_g = [](const vector_bigz& eval_pows) -> std::vector<G1> {
        // TODO parallelise by constructing res{eval_pows.size()}
        std::vector<G1> res{eval_pows.size()};
        assert(res.capacity() == 2 * N_);
        assert(res.size() == 2 * N_);
        #pragma omp parallel for schedule(static) num_threads(N_THREADS)
        for (size_t i = 0; i < res.size(); i++)
            res.at(i) = pow_(gen1, eval_pows[i]);
        assert(res.size() == 2 * N_);
        assert(res.capacity() == 2 * N_);
        return res;
    };
    // rounds half away from zero, e.g. 0.5 |-> 1 and -0.5 |-> -1
    auto round_vec = [](const vector_double& vec) -> vector_double {
        vector_double res(vec.size());
        for (size_t i = 0; i < vec.size(); i++)
            res[i] = mpf_round(vec[i]);
        return res;
    };
    auto map_to_q = [](vector_double& vec) -> vector_bigz {
        vector_bigz res(vec.size());
        for (size_t i = 0; i < vec.size(); i++) {
            mpf_class x = vec[i];
            bigz x1 = bigz(x);
            bigz val = x1 < 0 ? x1 + FIELD_MODULUS : x1;
            res.at(i) = val;
        }
        return res;
    };
    auto map_to_half_q = [](vector_bigz& vec) -> vector_double {
        vector_double res;
        for (bigz& x : vec) {
            mpf_class val = x > (FIELD_MODULUS / 2) ? x - FIELD_MODULUS : x;
            res.push_back(val);
        }
        return res;
    };

    Params pms;
    pms.print();

    DEBUG(std::cout << "*** START run_control_loop ***\n";)
    using matrix_i128 = array2d<bigz>;

    bigz q = pms.q;
    matrix_double A = pms.A;
    matrix_double B = pms.B;
    matrix_double C = pms.C;
    matrix_i128 F = pms.F;
    matrix_i128 G_bar = pms.G_bar;
    matrix_i128 H_bar = pms.H_bar;
    matrix_i128 R_bar = pms.R_bar;
    vector_double x_pl = pms.x_plant_init;
    vector_bigz x_cont = pms.x_cont_init_scaled;
    double rr = pms.r;
    double ss = pms.s;
    double L = pms.L;
    size_t iter_ = pms.iter_;
    DEBUG(iter_ = 3;)

    bigz from = 1;
    bigz to_inclusive = q - 1;
    // TODO Params should sample
    auto knowledge_exps = pms.sample_knowledge_exponents(from, to_inclusive);
    bigz alpha_0 = knowledge_exps.at(0);
    bigz alpha_1 = knowledge_exps.at(1);
    bigz gamma_0 = knowledge_exps.at(2);
    bigz gamma_1 = knowledge_exps.at(3);
    bigz rho_0 = knowledge_exps.at(4);
    bigz rho_1 = knowledge_exps.at(5);
    size_t m = H_bar.n_rows();
    size_t n = F.n_rows();
    // TODO Params should sample
    auto verification_vectors = pms.sample_verification_vectors(m, n, from, to_inclusive);
    auto r_0 = verification_vectors.at(0);
    auto r_1 = verification_vectors.at(1);
    auto s = verification_vectors.at(2);

    size_t N = N_;
    u32 d = pms.d;
    u32 power = pms.power;
    bigz v = pms.v;

    Encryptor enc(v, d, N, q);
    rgsw_mat F_ctx = enc.encrypt_rgsw_mat(F);
    rgsw_mat G_bar_ctx = enc.encrypt_rgsw_mat(G_bar);
    rgsw_mat R_bar_ctx = enc.encrypt_rgsw_mat(R_bar);
    rgsw_mat H_bar_ctx = enc.encrypt_rgsw_mat(H_bar);
    rlwe_vec x = enc.encrypt_rlwe_vec(x_cont);

    // TODO sample
    bigz eval_point = 42;
    vector_bigz eval_pows = eval_poly_pows(2 * N, eval_point, q);
    std::vector<G1> eval_pows_g = make_eval_pows_g(eval_pows);
    auto keys = compute_eval_and_veri_keys(
        F_ctx, G_bar_ctx, R_bar_ctx, H_bar_ctx,
        r_0, r_1, s, rho_0, rho_1, alpha_0, alpha_1, gamma_0, gamma_1,
        N, eval_pows, enc
    );
    auto ek = std::get<0>(keys);
    veri_key vk = std::get<1>(keys);
    check_key ck = std::get<2>(keys);

    // convert rgsw mats to ntt/eval form
    std::cout << "start eval\n";
    F_ctx.to_eval_form();
    G_bar_ctx.to_eval_form();
    R_bar_ctx.to_eval_form();
    H_bar_ctx.to_eval_form();
    std::cout << "end eval\n";

    // create initial proof
    auto x_hash = x.get_hash(eval_pows); // XXX
    // TODO move to veri_vec
    bigz rx_0 = r_1.dot_prod(x_hash); // XXX
    GT g1 = pow_t(genT, rx_0);
    Proof old_proof {};
    old_proof.grx_ = g1;

    TIMING(timing.iter_ = 1;)
    TIMING(print_times_and_counts(timing);)
    TIMING(timing = {};)
    TIMING(timing.iter_ = iter_;)
    TIMING(auto start = std::chrono::high_resolution_clock::now();)
    for (size_t k = 0; k < iter_; k++) {
        DEBUG(std::cout << "\n*** START k: " << k << ", run_control_loop ***\n";)
        TIMING(auto start_plant = std::chrono::high_resolution_clock::now();)

        /*  ### Plant: Compute output ### */
        vector_double y_out = mat_vec_mult(C, x_pl);
        vars.y.push_back(y_out);
        vector_double y_out_scaled = scalar_vec_mult(rr, y_out);
        vector_double y_out_rounded = round_vec(y_out_scaled);
        y_out_scaled = scalar_vec_mult(L, y_out_rounded);
        y_out_rounded = round_vec(y_out_scaled);
        vector_bigz y_out_q = map_to_q(y_out_rounded);
        rlwe_vec y = enc.encrypt_rlwe_vec(y_out_q); // P -> C

        #ifdef TIMING_ON
            auto end_plant = std::chrono::high_resolution_clock::now();
            timing.plant += end_plant - start_plant;
            auto start_controller = std::chrono::high_resolution_clock::now();
        #endif

        /* ### Controller: Compute output ### */
        // TODO could create a copy for passing to compute_proof() to save
        // re-creation of x_d
        rlwe_decomp_vec x_d = x.decompose(v, d, power);
        x_d.to_eval_form();
        rlwe_vec u_conv = H_bar_ctx.convolve(x_d); // C -> P
        u_conv.to_coeff_form();

        #ifdef TIMING_ON
            auto end_controller = std::chrono::high_resolution_clock::now();
            timing.controller += end_controller - start_controller;
            start_plant = std::chrono::high_resolution_clock::now();
        #endif

        /* ### Plant: Update state (and re-encrypt u) ### */
        rlwe_vec u = u_conv.mod_cyclo(N);
        vector_bigz u_ptx = enc.decrypt_rlwe_vec(u);
        vector_double u_ptx_q2 = map_to_half_q(u_ptx);
        vector_double u_in = scalar_vec_mult(1.0 / (rr * ss * ss * L), u_ptx_q2);
        // vars.u.push_back(u_in); // XXX
        vector_double u_in_scaled = round_vec(scalar_vec_mult(rr, u_in));
        u_in_scaled = round_vec(scalar_vec_mult(L, u_in_scaled));
        vector_bigz u_in_q = map_to_q(u_in_scaled);
        rlwe_vec u_re = enc.encrypt_rlwe_vec(u_in_q); // XXX P->C: Plant re-encrypts u

        x_pl = mat_vec_mult(A, x_pl);
        vector_double B_u = mat_vec_mult(B, u_in);
        for (size_t i = 0; i < B_u.size(); i++)
            x_pl.at(i) += B_u.at(i);
        // vars.x_plant.push_back(x_pl); // XXX

        #ifdef TIMING_ON
            end_plant = std::chrono::high_resolution_clock::now();
            timing.plant += end_plant - start_plant;
            start_controller = std::chrono::high_resolution_clock::now();
        #endif

        /* ###  Controller: Update state ### */
        // rlwe_vec F_x = F_ctx.convolve(x_d);
        // rlwe_decomp_vec y_d = y.decompose(v, d, power);
        // y_d.to_eval_form();
        // rlwe_vec G_y = G_bar_ctx.convolve(y_d);
        // rlwe_decomp_vec u_re_d = u_re.decompose(v, d, power);
        // u_re_d.to_eval_form();
        // rlwe_vec R_u = R_bar_ctx.convolve(u_re_d);
        // rlwe_vec x_conv = F_x + G_y + R_u;

        rlwe_vec x_conv = F_ctx.convolve(x_d) \
            + G_bar_ctx.convolve(y.decompose(v, d, power).to_eval_form()) \
            + R_bar_ctx.convolve(u_re.decompose(v, d, power).to_eval_form());
        x_conv.to_coeff_form();
        rlwe_vec x_old = x;
        x = x_conv.mod_cyclo(N); // FIXME remove N

        // TODO create x_hi
        rlwe_vec x_hi{x_conv.n_rlwes(), x_conv.n_polys(), N_};
        for (size_t i = 0; i < x_hi.n_rlwes(); i++) {
            for (size_t j = 0; j < x_hi.n_polys(); j++) {
                #pragma omp parallel for schedule(static) num_threads(N_THREADS)
                for (size_t k = 0; k < N_; k++)
                    x_hi[i][j][k] = x_conv[i][j][k + N_];
            }
        }
        // std::cout << "x_hi[3][1][N_ - 1]: " << x_hi[3][1][N_ - 1] << "\n";

        // // Decrypt and capture controller state
        // vector_bigz x_ptx = enc.decrypt_rlwe_vec(x);
        // vector_double x_ptx_q2 = map_to_half_q(x_ptx);
        // vector_double x_ptx_scaled = scalar_vec_mult(1 / (rr * ss * L), x_ptx_q2);
        // vars.x_cont.push_back(x_ptx_scaled);

        #ifdef TIMING_ON
            auto start_proof = std::chrono::high_resolution_clock::now();
        #endif
        // std::cout << "before prove\n";
        /* ### Controller: Prove ### */
        const eval_key& ek_i = (k % 2 == 0) ? std::get<0>(ek) : std::get<1>(ek);
        // TODO x_conv->x_hi
        // Proof proof = compute_proof(ek_i, x, x_conv, x_old, v, d, power, eval_pows_g); // C -> P
        Proof proof = compute_proof(ek_i, x, x_hi, x_old, v, d, power, eval_pows_g); // C -> P
        // std::cout << "after prove\n";

        #ifdef TIMING_ON
            auto end_proof = std::chrono::high_resolution_clock::now();
            timing.proof += end_proof - start_proof;
            end_controller = std::chrono::high_resolution_clock::now();
            timing.controller += end_controller - start_controller;
            start_plant = std::chrono::high_resolution_clock::now();
            auto start_verify = std::chrono::high_resolution_clock::now();
        #endif

        /* ### Plant: Verify ### */
        // std::cout << "before veri\n";
        verify_with_lin_and_dyn_checks(vk, proof, old_proof, k, y, u_conv, u_re, v, d, power, eval_pows);
        // std::cout << "after veri\n";
        old_proof = proof;

        #ifdef TIMING_ON
            auto end_verify = std::chrono::high_resolution_clock::now();
            timing.verify += end_verify - start_verify;
            end_plant = std::chrono::high_resolution_clock::now();
            timing.plant += end_plant - start_plant;
        #endif

        DEBUG(std::cout << "\n\n*** END k: " << k << ", run_control_loop ***\n\n";)
    }
    TIMING(auto end = std::chrono::high_resolution_clock::now();)
    TIMING(timing.loop = end - start;)
}
void run_control_loop_unencrypted(control_law_vars& vars) {
    Params pms;
    // DEBUG(pms.print();)
    DEBUG(std::cout << "*** START run_control_loop_unencrypted ***\n";)
    // using matrix_i128 = array2d<bigz>;

    matrix_double A = pms.A;
    matrix_double B = pms.B;
    matrix_double C = pms.C;
    matrix_double F = pms.get_F_();
    matrix_double G = pms.get_G();
    matrix_double H = pms.get_H();
    matrix_double R = pms.get_R();
    vector_double x_plant_init = pms.x_plant_init;
    vector_double x_cont_init = pms.get_x_cont_init();
    size_t iter_ = pms.iter_;
    DEBUG(iter_ = 3;)

    DEBUG(std::cout << "A:\n";)
    DEBUG(for (auto& x : A) print_vector_double(x);)
    DEBUG(std::cout << "B:\n";)
    DEBUG(for (auto& x : B) print_vector_double(x);)
    DEBUG(std::cout << "C:\n";)
    DEBUG(for (auto& x : C) print_vector_double(x);)
    DEBUG(std::cout << "F:\n";)
    DEBUG(for (auto& x : F) print_vector_double(x);)
    DEBUG(std::cout << "G:\n";)
    DEBUG(for (auto& x : G) print_vector_double(x);)
    DEBUG(std::cout << "H:\n";)
    DEBUG(for (auto& x : H) print_vector_double(x);)

    DEBUG(std::cout << "x_plant_init:\n";)
    DEBUG(print_vector_double(x_plant_init);)
    DEBUG(std::cout << "x_cont_init:\n";)
    DEBUG(print_vector_double(x_cont_init);)

    // Initialise Plant and Controller states
    vector_double x_plant(x_plant_init.size());
    vector_double x_cont(x_cont_init.size());
    assert(x_plant_init.size() == x_cont_init.size());
    for (size_t i = 0; i < x_plant_init.size(); i++) {
        x_plant.at(i) = x_plant_init.at(i);
        x_cont.at(i) = x_cont_init.at(i);
    }
    DEBUG(std::cout << "x_plant:\n";)
    DEBUG(print_vector_double(x_plant);)
    DEBUG(std::cout << "x_cont:\n";)
    DEBUG(print_vector_double(x_cont);)


    for (size_t k = 0; k < iter_; k++) {
        DEBUG(std::cout << "\n*** START k: " << k << ", run_control_loop_unenc ***\n";)

        // Plant outputs y = Cx
        vector_double y = mat_vec_mult(C, x_plant);
        vars.y.push_back(y);
        DEBUG(std::cout << "y:\n";)
        DEBUG(print_vector_double(y);)

        // Controller outputs u = Hx
        vector_double u = mat_vec_mult(H, x_cont);
        vars.u.push_back(u);
        DEBUG(std::cout << "u:\n";)
        DEBUG(print_vector_double(u);)

        // Plant updates state x' = Ax + Bu
        vector_double Ax = mat_vec_mult(A, x_plant);
        vector_double Bu = mat_vec_mult(B, u);
        assert(Ax.size() == Bu.size());
        assert(x_plant.size() == Bu.size());
        for (size_t i = 0; i < x_plant.size(); i++)
            x_plant.at(i) = Ax.at(i) + Bu.at(i);
        vars.x_plant.push_back(x_plant);
        // DEBUG(std::cout << "x_plant:\n";)
        // DEBUG(print_vector_double(x_plant);)

        // Controller updates state x' = Fx + Gy
        vector_double Fx = mat_vec_mult(F, x_cont);
        vector_double Gy = mat_vec_mult(G, y);
        vector_double Ru = mat_vec_mult(R, u);
        assert(Fx.size() == Gy.size());
        assert(Ru.size() == Gy.size());
        assert(x_cont.size() == Gy.size());
        for (size_t i = 0; i < x_cont.size(); i++)
            x_cont.at(i) = Fx.at(i) + Gy.at(i) + Ru.at(i);
        vars.x_cont.push_back(x_cont);
        // DEBUG(std::cout << "x_cont:\n";)
        DEBUG(print_vector_double(x_cont);)

        DEBUG(std::cout << "\n\n*** END k: " << k << ", run_control_loop_unenc ***\n\n";)
    }
}

void print_vars_diff(control_law_vars& vars, control_law_vars& vars_unenc) {
    std::cout << "##### VARS_UNENC #####\n";
    std::cout << "y:\n";
    for (auto x : vars_unenc.y)
        print_vector_double(x);
    std::cout << "u:\n";
    for (auto x : vars_unenc.u)
        print_vector_double(x);
    std::cout << "x_plant:\n";
    for (auto x : vars_unenc.x_plant)
        print_vector_double(x);
    std::cout << "x_cont:\n";
    for (auto x : vars_unenc.x_cont)
        print_vector_double(x);
    std::cout << "##### VARS #####\n";
    std::cout << "y:\n";
    for (auto x : vars.y)
        print_vector_double(x);
    std::cout << "u:\n";
    for (auto x : vars.u)
        print_vector_double(x);
    std::cout << "x_plant:\n";
    for (auto x : vars.x_plant)
        print_vector_double(x);
    std::cout << "x_cont:\n";
    for (auto x : vars.x_cont)
        print_vector_double(x);
    // control_law_vars diffs;
    for (size_t i = 0; i < vars.y.size(); i++)
        for (size_t j = 0; j < vars.y.at(0).size(); j++)
            vars_unenc.y.at(i).at(j) -= vars.y.at(i).at(j);
    for (size_t i = 0; i < vars.u.size(); i++)
        for (size_t j = 0; j < vars.u.at(0).size(); j++)
            vars_unenc.u.at(i).at(j) -= vars.u.at(i).at(j);
    for (size_t i = 0; i < vars.x_plant.size(); i++)
        for (size_t j = 0; j < vars.x_plant.at(0).size(); j++)
            vars_unenc.x_plant.at(i).at(j) -= vars.x_plant.at(i).at(j);
    for (size_t i = 0; i < vars.x_cont.size(); i++)
        for (size_t j = 0; j < vars.x_cont.at(0).size(); j++)
            vars_unenc.x_cont.at(i).at(j) -= vars.x_cont.at(i).at(j);


    std::cout << "##### DIFFS #####\n";
    std::cout << "y:\n";
    for (auto x : vars_unenc.y)
        print_vector_double(x);
    std::cout << "u:\n";
    for (auto x : vars_unenc.u)
        print_vector_double(x);
    std::cout << "x_plant:\n";
    for (auto x : vars_unenc.x_plant)
        print_vector_double(x);
    std::cout << "x_cont:\n";
    for (auto x : vars_unenc.x_cont)
        print_vector_double(x);


}

int main() {
    // sleep(20);
    control_law_vars vars;
    control_law_vars vars_unenc;
    // omp_set_nested(1);
    TIMING(auto start = std::chrono::high_resolution_clock::now();)

    initPairing(BN_SNARK1);

    int gen_seed = 42;
    hashAndMapToG1(gen1, std::string("P_") + std::to_string(gen_seed));
    hashAndMapToG2(gen2, std::string("P_") + std::to_string(gen_seed));
    pairing(genT, gen1, gen2);

    // std::cout << "sizeof (long int)=" << sizeof (long int) << "\n";
    run_control_loop(vars, timing);

    #ifdef TIMING_ON
        auto end = std::chrono::high_resolution_clock::now();
        timing.total = end - start;
        print_times_and_counts(timing);
    #endif

    // std::cout << std::fixed << std::setprecision(2);
    // run_control_loop_unencrypted(vars_unenc);
    // print_vars_diff(vars, vars_unenc);
    return 0;
}
