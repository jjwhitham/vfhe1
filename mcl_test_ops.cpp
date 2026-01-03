#include <gmpxx.h>
#include "serialisation.hpp"
#include "/Users/jw/Projects/mcl/include/mcl/bn.hpp"
#include <iostream>
#include <vector>

#ifndef N_THREADS
    #define N_THREADS 1
#endif

typedef mpz_class bigz;
using namespace mcl::bn;

bigz FIELD_MODULUS("\
21888242871839275222246405745257275088548364400416034343698204186575808495617\
");

G1 pow_(const G1& base, const bigz& power) {
    assert(power >= 0);
    // assert(power < FIELD_MODULUS);
    // static Fr power1;
    Fr power1;
    power1.clear();
    if (power > FIELD_MODULUS) {
        std::cout << "pow_: power > FIELD_MODULUS\n";
        mpz_to_Fr(power1, power % FIELD_MODULUS);
    }
    else
        mpz_to_Fr(power1, power);
    return base * power1;
}

G2 pow2_(const G2& base, const bigz& power) {
    assert(power >= 0);
    // assert(power < FIELD_MODULUS);
    // static Fr power1;
    Fr power1;
    power1.clear();
    if (power > FIELD_MODULUS) {
        std::cout << "pow_2: power > FIELD_MODULUS\n";
        mpz_to_Fr(power1, power % FIELD_MODULUS);
    }
    else
        mpz_to_Fr(power1, power);
    return base * power1;
}

GT pow_t(const GT& base, const bigz& power) {
    assert(power >= 0);
    // assert(power < FIELD_MODULUS);
    // static Fr power1;
    Fr power1;
    power1.clear();
    if (power > FIELD_MODULUS) {
        std::cout << "pow_t: power > FIELD_MODULUS\n";
        mpz_to_Fr(power1, power % FIELD_MODULUS);
    }
    else
        mpz_to_Fr(power1, power);
    GT ret;
    // pairing(ret, Generator, Gen2); // TODO move out into main() and call it GenT
    GT::pow(ret, base, power1);
    return ret;
}

void test_powers_pairing(G1 g1, G2 g2) {
    // Tests the pow and pairing funcs, asserting:
    // e(g1, g2) == e(g1^a, g2^b) == e(g1, g2)^ab

    // Create scalars
    bigz a{42};
    bigz b{666};
    bigz ab = a * b;

    // do pairing 1
    GT p1;
    // p1.clear();
    pairing(p1, pow_(g1, a), pow2_(g2, b));

    // do pairing 2
    GT p2;
    // p2.clear();
    pairing(p2, g1, pow2_(g2, ab));

    // do pairing 3
    GT p3;
    // p3.clear();
    pairing(p3, pow_(g1, ab), g2);

    // do pairing 4
    GT p4;
    // p4.clear();
    pairing(p4, g1, g2);
    p4 = pow_t(p4, ab);

    // std::cout << "p1: " << p1 << "\n\n";
    // std::cout << "p2: " << p2 << "\n\n";
    // std::cout << "p3: " << p3 << "\n\n";
    // std::cout << "p4: " << p4 << "\n\n";

    assert(p1 == p2);
    assert(p1 == p3);
    assert(p1 == p4);
}

G1 msm(std::vector<bigz>& self, std::vector<G1>& eval_pows_g) {
    // for each of the pols in rlwe_decomp
    // msm creates a vector of polys converted to Fr's
    size_t N_ = 4096;
    assert(self.size() <= 2 * N_);
    std::vector<G1> eval_pows_g1;
    eval_pows_g1.reserve(self.size());
    for (size_t i = 0; i < self.size(); i++)
        eval_pows_g1.push_back(eval_pows_g.at(i));

    std::vector<Fr> scalars(self.size());
    // #pragma omp parallel for schedule(static) num_threads(N_THREADS)
    for (size_t i = 0; i < self.size(); i++)
        scalars.at(i).clear();

    // #pragma omp parallel for schedule(static) num_threads(N_THREADS)
    for (size_t i = 0; i < self.size(); i++) {
        mpz_to_Fr(scalars[i], self.at(i));
    }

    G1 res;
    G1::mulVecMT(res, eval_pows_g1.data(), scalars.data(), self.size(), N_THREADS);
    return res;
}

bigz hash(const std::vector<bigz>& poly, const bigz& t) {
    bigz t_pow{1};
    bigz hash_val{0};
    for (size_t i = 0; i < poly.size(); i++) {
        hash_val += (poly.at(i) * t_pow) % FIELD_MODULUS;
        t_pow = (t_pow * t) % FIELD_MODULUS;
    }
    hash_val %= FIELD_MODULUS;
    return hash_val;
}

void test_msm(G1 g1, G2 g2) {
    // create a vector of power's of t
    size_t TWO_N = 8192;
    std::vector<G1> powers_of_t_;
    powers_of_t_.reserve(TWO_N);
    assert(powers_of_t_.capacity() == TWO_N);
    bigz t_pow{1};
    int t = 42;
    for (size_t i = 0; i < powers_of_t_.capacity(); i++) {
        G1 g1_pow = pow_(g1, t_pow);
        powers_of_t_.push_back(g1_pow);
        t_pow = (t_pow * t) % FIELD_MODULUS;
    }
    assert(powers_of_t_.size() == 8192);

    // random polys
    bigz a_val{42};
    std::vector<bigz> a;
    a.reserve(TWO_N);
    bigz b_val{666};
    std::vector<bigz> b;
    b.reserve(TWO_N);
    for (size_t i = 0; i < TWO_N; i++) {
        a.push_back(a_val);
        b.push_back(b_val);
        a_val++;
        b_val++;
    }
    std::vector<G1> g1_a;
    g1_a.reserve(TWO_N);
    for (size_t i = 0; i < TWO_N; i++) {
        G1 g1_val = pow_(g1, a.at(i));
        g1_a.push_back(g1_val);
    }
    // hash H(ab) = H(a)H(b)
    bigz hash_a = hash(a, t);
    bigz hash_b = hash(b, t);
    bigz hash_ab = hash_a * hash_b;
    // pairing(g1, g2)^H(ab)
    GT p1;
    p1.clear();
    pairing(p1, g1, g2);
    p1 = pow_t(p1, hash_ab);

    // msm(a)
    G1 msm_a = msm(a, powers_of_t_);
    // g2^H(b)
    G2 g2_hash_b = pow2_(g2, hash_b);
    // pairing(msm(a), g2^H(b))
    GT p2;
    p2.clear();
    pairing(p2, msm_a, g2_hash_b);
    std::cout << "p1: " << p1 << "\n";
    std::cout << "p2: " << p2 << "\n";
    assert(p1 == p2);
}

void test_gt(G1 g1, G2 g2) {
    GT res{1};
    std::cout << "\n\nres: " << res << "\n";
    res.clear();
    std::cout << "res: " << res << "\n";

    Fr num = 1;
    GT genT;
    pairing(genT, g1, g2);
    std::cout << "\ngenT: " << genT << "\n";
    GT res1 = pow_t(res, 0);
    std::cout << "\n\nres1: " << res1 << "\n";
}

int main() {
    initPairing(BN_SNARK1);

    // create g1 gen
    int gen_seed = 42;

    G1 g1;
    hashAndMapToG1(g1, std::string("P_") + std::to_string(gen_seed));

    // create g2 gen
    G2 g2;
    hashAndMapToG2(g2, std::string("P_") + std::to_string(gen_seed));

    test_powers_pairing(g1, g2);
    test_msm(g1, g2);
    test_gt(g1, g2);

    return 0;
}