#include "vfhe.h"
#include <chrono>

void init(rgsw_mat& F, rlwe_decomp_vec& x, veri_vec_scalar& r) {
    std::random_device rd;  // Non-deterministic seed
    std::mt19937 gen(rd()); // Mersenne Twister engine
    std::uniform_int_distribution<i128> distrib(0, FIELD_MODULUS - 1); // Range: 0 to 100
    i128 random_number = distrib(gen);
    i128 counter = random_number;
    for (size_t i = 0; i < F.n_rows(); ++i) {
        for (size_t j = 0; j < F.n_cols(); ++j) {
            rgsw& rg = F.get_rgsw(i, j);
            for (rlwe& rl : rg) {
                for (poly& p : rl) {
                    for (size_t m = 0; m < p.size(); ++m) {
                        p.set(m, counter++ % FIELD_MODULUS); // just an example initialization
                    }
                }
            }
        }
    }
    counter = distrib(gen);
    for (size_t i = 0; i < x.size(); ++i) {
        rlwe_decomp& rd = x.get(i);
        for (size_t j = 0; j < rd.size(); ++j) {
            poly& p = rd.get(j);
            for (size_t m = 0; m < p.size(); ++m) {
                p.set(m, counter++ % FIELD_MODULUS);
            }
        }
    }
    counter = distrib(gen);
    for (size_t i = 0; i < r.size(); ++i) {
        r.set(i, counter++ % FIELD_MODULUS);
    }
}

void test() {
    rlwe_vec u(2, 2, 2);
    for (size_t i = 0; i < u.size(); i++) {
        auto& x = u.get(i);
        for (size_t j = 0; j < x.size(); j++) {
            auto& y = x.get(j);
            for (size_t k = 0; k < y.size(); k++) {
                y.set(k, k);
            }
        }
    }
    std::cout << "u:\n";
    u.print();
    veri_vec_scalar r(2);
    for (size_t i = 0; i < r.size(); i++) {
        r.set(i, i + 1);
    }
    std::cout << "r:\n";
    r.print();
    std::cout << "\n";
    auto ru = r * u;

    std::cout << "ru:\n";
    ru.print();
    i128 expected[] = {0, 3, 6, 9};
    for (auto& x : ru) {
        for (size_t i = 0; i < x.size(); i++) {
            assert(x.get(i) == expected[i]);
        }
    }
    i128 expected1[] = {1, 18, 2, 13};
    rlwe gru = ru.pow();
    std::cout << "gru:\n";
    gru.print();
    for (auto& x : gru) {
        for (size_t i = 0; i < x.size(); i++) {
            assert(x.get(i) == expected1[i]);
        }
    }
    veri_vec_scalar gr = r.pow();
    std::cout << "gr:\n";
    gr.print();
    std::cout << "\n";
    i128 expected2[] = {4, 16};
    for (size_t i = 0; i < r.size(); i++) {
        assert(gr.get(i) == expected2[i]);
    }
    rlwe_vec gu = u.pow();
    std::cout << "gu:\n";
    gu.print();
    i128 expected3[] = {1, 4, 16, 18};
    for (auto& rlwe_el : gu)
        for (auto& poly_el : rlwe_el)
            for (size_t i = 0; i < poly_el.size(); i++) {
                assert(poly_el.get(i) == expected3[i]);
            }

    std::cout << "gru1:\n";
    rlwe gru1 = gr.pow(u);
    gru1.print();

    poly p(2);
    p.set(0, 1);
    p.set(1, 2);
    rlwe_decomp rd(2, 2);
    for (size_t i = 0; i < rd.size(); i++)
        rd.set(i, p); // set all decompositions to the same poly
    std::cout << "rdv:\n";
    rd.print();
    rd.get(0).set(0, 3); // change first coeff of first poly
    std::cout << "rdv after change:\n";
    rd.print();
}

void test_rlwe_decomp() {
    size_t n_coeffs = 2;
    size_t n_polys = 2;
    size_t n_rlwes = 2;
    rlwe_decomp rd(n_polys, n_coeffs);
    for (size_t i = 0; i < rd.size(); i++) {
        poly p(n_coeffs);
        for (size_t j = 0; j < p.size(); j++) {
            p.set(j, j + 1);
        }
        rd.set(i, p);
    }
    std::cout << "rlwe_decomp:\n";
    rd.print();
    rgsw rg(n_rlwes, n_polys, n_coeffs);
    for (size_t i = 0; i < rg.size(); i++) {
        rlwe r(n_polys, n_coeffs);
        for (size_t j = 0; j < r.size(); j++) {
            poly p(n_coeffs);
            for (size_t k = 0; k < p.size(); k++) {
                p.set(k, k + 1);
            }
            r.set(j, p);
        }
        rg.set(i, r);
    }
    std::cout << "rgsw:\n";
    rg.print();

    rlwe rl = rg * rd; // rgsw * rlwe_decomp
    std::cout << "rlwe from rgsw * rlwe_decomp:\n";
    rl.print();

    i128 r = 2;
    rlwe rrl = rl * r; // scalar * rlwe
    std::cout << "rlwe from scalar * rlwe:\n";
    rrl.print();


    rlwe grrl = rrl.pow();
    std::cout << "grrl";
    grrl.print();
    i128 gr = 16; // g^r = 4^2 % 23 = 16
    rlwe grl = rl.pow(gr); // g^r.pow(rlwe)
    std::cout << "grl:\n";
    grl.print();
}

void test_rlwe_decomp_vec() {
    size_t n_coeffs = 2;
    size_t n_polys = 2;
    size_t n_rlwes = 2;
    size_t n_rlwe_decomps = 2;
    rlwe_decomp_vec x(n_rlwe_decomps, n_polys, n_coeffs);
    i128 counter = 0;
    for (size_t v = 0; v < x.size(); v++) {
        rlwe_decomp& rd = x.get(v);
        for (size_t i = 0; i < rd.size(); i++) {
            poly& p = rd.get(i);
            for (size_t j = 0; j < p.size(); j++) {
                p.set(j, counter++ % FIELD_MODULUS); // initialize coefficients to j + 1
            }
            rd.set(i, p);
        }
        x.set(v, rd);
    }
    std::cout << "x:\n";
    x.print();
    size_t rows = 2;
    size_t cols = 2;
    counter = 0;
    rgsw_mat F(rows, cols, n_rlwes, n_polys, n_coeffs);
    for (size_t m = 0; m < F.n_rows(); m++) {
        for (size_t n = 0; n < F.n_cols(); n++) {
            rgsw& rg = F.get(m, n);
            for (size_t i = 0; i < rg.size(); i++) {
                rlwe& r = rg.get(i);
                for (size_t j = 0; j < r.size(); j++) {
                    poly& p = r.get_poly(j);
                    for (size_t k = 0; k < p.size(); k++) {
                        p.set(k, (counter++) % FIELD_MODULUS);
                    }
                }
            }
        }
    }
    std::cout << "F:\n";
    for (size_t i = 0; i < F.n_rows(); i++) {
        for (size_t j = 0; j < F.n_cols(); j++) {
            F.get(i, j).print();
            std::cout << "\n";
        }
    }

    std::cout << "Fx:\n";
    rlwe_vec Fx = F * x; // rgsw_mat * rlwe_decomp_vec
    Fx.print();

    veri_vec_scalar r(n_rlwes);
    for (size_t i = 0; i < r.size(); i++) {
        r.set(i, i + 1);
    }
    std::cout << "r:\n";
    r.print();

    rlwe rFx = r * Fx; // veri_vec_scalar * rlwe_vec
    std::cout << "rFx:\n";
    rFx.print();

    rgsw_vec rF = r * F; // veri_vec_scalar * rgsw_mat
    std::cout << "rF:\n";
    rF.print();


    rgsw_vec grF = rF.pow(); // g^r
    std::cout << "grF:\n";
    grF.print();
    rlwe grFx = grF.pow(x); // g^r.pow(Fx)
    std::cout << "grFx:\n";
    grFx.print();

    rlwe rFx1 = rF * x; // rgsw_vec * rlwe_decomp_vec
    std::cout << "rFx1:\n";
    rFx1.print();

    rlwe grFx1 = rFx1.pow(); // g^rFx1
    std::cout << "grFx1:\n";
    grFx1.print();

    for (size_t i = 0; i < grFx.size(); i++) {
        poly& p1 = grFx.get_poly(i);
        poly& p2 = grFx1.get_poly(i);
        for (size_t j = 0; j < p1.size(); j++) {
            assert(p1.get(j) == p2.get(j)); // check if grFx and grFx1 are equal
        }
    }
}

void test_full() {
    size_t n_rlwes = 7;
    size_t n_polys = 2; // XXX must be 2
    size_t n_coeffs = 4096;
    size_t rows = 5;
    size_t cols = 12;
    rgsw_mat F(rows, cols, n_rlwes, n_polys, n_coeffs);
    rlwe_decomp_vec x(cols, n_rlwes, n_coeffs);
    veri_vec_scalar r(rows);
    init(F, x, r);
    // for (size_t i = 0; i < F.n_rows(); i++) {
    //     for (size_t j = 0; j < F.n_cols(); j++) {
    //         F.get(i, j).print();
    //         std::cout << "\n";
    //     }
    // }
    // x.print();
    // r.print();
    // std::cout << "\n";

    // Compute Fx = F * x
    rlwe_vec Fx = F * x;

    // Compute rF = r * F
    rgsw_vec rF = r * F;

    // Compute grF = g^{rF}
    rgsw_vec grF = rF.pow();

    auto start = std::chrono::high_resolution_clock::now();
    // Compute grFx = (g^{rF})^x = g^{rF * x}
    rlwe grFx = grF.pow(x);
    auto end = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double, std::nano> elapsed = end - start;
    std::cout << "Elapsed time: " << elapsed.count() << " ns\n";
    // init(F, x, r); // reinitialize F, x, r

    // Compute rFx = rF * x
    rlwe rFx = rF * x;
    // rlwe rFx1 = r * Fx;

    // Compute grFx2 = g^{rFx}
    rlwe grFx2 = rFx.pow();

    // init(F, x, r); // reinitialize F, x, r

    // Compute gr = g^r (elementwise)
    veri_vec_scalar gr = r.pow();

    // Compute grFx3 = (g^r)^{Fx} = g^{r * Fx}
    rlwe grFx3 = gr.pow(Fx);

    // std::cout << "grFx:\n";
    // for (const auto& poly_el : grFx) {
    //     for (const auto& coeff : poly_el)
    //         std::cout << i128str(coeff) << ", ";
    //     std::cout << "\n";
    // }
    // std::cout << "grFx2:\n";
    // for (const auto& poly_el : grFx2) {
    //     for (const auto& coeff : poly_el)
    //         std::cout << i128str(coeff) << ", ";
    //     std::cout << "\n";
    // }
    // std::cout << "grFx3:\n";
    // for (const auto& poly_el : grFx3) {
    //     for (const auto& coeff : poly_el)
    //         std::cout << i128str(coeff) << ", ";
    //     std::cout << "\n";
    // }
    // // std::cout << "r:\n";
    // // r.print();
    // // std::cout << "F:\n";
    // // for (size_t i = 0; i < F.n_rows(); i++) {
    // //     for (size_t j = 0; j < F.n_cols(); j++) {
    // //         F.get(i, j).print();
    // //         std::cout << "\n";
    // //     }
    // // }
    // // std::cout << "rF:\n";
    // // rF.print();

    // assert coeffs are equal for grFx, grFx2, grFx3
    for (size_t i = 0; i < grFx.size(); i++) {
        poly& p1 = grFx.get_poly(i);
        poly& p2 = grFx2.get_poly(i);
        poly& p3 = grFx3.get_poly(i);
        for (size_t j = 0; j < p1.size(); j++) {
            assert(p1.get(j) == p2.get(j));
            assert(p1.get(j) == p3.get(j));
        }
    }

}

int main() {
    test();
    test_rlwe_decomp();
    test_rlwe_decomp_vec();
    test_full();
    // omp_set_nested(1);
    // omp_set_num_threads(1);
    // std::cout << "omp num threads: " << omp_get_max_threads() << std::endl;
    // std::cout << " omp num threads: " << omp_get_num_threads() << std::endl;
    test_full();

    return 0;
}
