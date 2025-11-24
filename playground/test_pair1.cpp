// bench_bn254.cpp
#include <pbc/pbc.h>
// #include <gmp.h>
#include <gmpxx.h>
#include <chrono>
#include <iostream>
#include <vector>
#include <string>
#include <memory>

using mpz = mpz_class;

// RAII wrapper for PBC element_t (stack-allocated array-of-1 idiom)
class PBCElement {
public:
    enum Type { G1, G2, GT };
private:
    element_t e;         // array-of-1 storage used by PBC
    bool inited = false;
public:
    PBCElement(pairing_t pairing, Type t) {
        switch (t) {
            case G1: element_init_G1(e, pairing); break;
            case G2: element_init_G2(e, pairing); break;
            case GT: element_init_GT(e, pairing); break;
        }
        inited = true;
    }
    ~PBCElement() { if (inited) element_clear(e); }

    // non-copyable, non-movable (safe and simple)
    PBCElement(const PBCElement&) = delete;
    PBCElement& operator=(const PBCElement&) = delete;
    PBCElement(PBCElement&&) = delete;
    PBCElement& operator=(PBCElement&&) = delete;

    // access underlying pointer expected by PBC API
    element_ptr ptr() { return e; }
    const element_ptr ptr() const { return const_cast<element_t&>(e); }
};

// Precomputed BN254-like parameter string (example). Replace with an exact trusted BN254 param if available.
const char *bn254_param = R"param(
type f
q 205523667896953300194896352429254920972540065223
r 205523667896953300194895899082072403858390252929
b 40218105156867728698573668525883168222119515413
beta 115334401956802802075595682801335644058796914268
alpha0 191079354656274778837764015557338301375963168470
alpha1 71445317903696340296199556072836940741717506375
)param";

/*
struct f_param_s {
    mpz_t q; // Curve defined over F_q.
    mpz_t r; // The order of the curve.
    mpz_t b; // E: y^2 = x^3 + b
    mpz_t beta; //beta is a quadratic nonresidue in Fq
        //we use F_q^2 = F_q[sqrt(beta)]
    mpz_t alpha0, alpha1;
        //the polynomial x^6 + alpha0 + alpha1 sqrt(beta)
        //is irreducible over F_q^2[x], so
        //we can extend F_q^2 to F_q^12 using the
        //sixth root of -(alpha0 + alpha1 sqrt(beta))
};
*/


int main() {
  pairing_t pairing;
  if (pairing_init_set_str(pairing, bn254_param) != 0) {
    std::cerr << "pairing_init_set_str failed\n";
    return 1;
  }

  // Initialize elements
//   element_t g1, g2, tG; // g1 in G1, g2 in G2, tG in GT
//   element_init_G1(g1, pairing);
//   element_init_G2(g2, pairing);
//   element_init_GT(tG, pairing);
  PBCElement g1(pairing, PBCElement::G1);
  PBCElement g2(pairing, PBCElement::G2);
  PBCElement tG(pairing, PBCElement::GT);

  // Set generators: use random for demonstration (or hash-to-point)
  element_random(g1.ptr());
  element_random(g2.ptr());

  // Prepare scalar mpz_t
  mpz_t scalar;
  mpz_init(scalar);
  mpz_t tmp;
  mpz_init(tmp);

  // Pre-generate 1M random small scalars (mod r) for exponentiations to avoid RNG cost in loop.
  const size_t N1 = 1000; // 1e6 for G1
  const size_t N2 = 1000; // 1e6 for G2
  const size_t NP = 100;    // 1e3 pairings

  std::vector<mpz> scalars1;
  std::vector<mpz> scalars2;
  scalars1.reserve(N1);
  scalars2.reserve(N2);

  // Get group order r as mpz_t
  mpz_t order;
  mpz_init(order);
  // pairing->r is internal; use element order via element_order_mpz if available.
  // Fallback: use mpz_set_str from param r value
  mpz_set_str(order, "21888242871839275222246405745257275088548364400416034343698204186575808495617", 10);

  // RNG: generate scalars by randomizing via mpz_urandomm with a gmp_randstate_t
  gmp_randstate_t state;
  gmp_randinit_default(state);
  gmp_randseed_ui(state, (unsigned long)time(nullptr));

  auto gen_scalars = [&order, &state](std::vector<mpz> &vec, size_t n) { // , gmp_randstate_t& state) {
    for (size_t i = 0; i < n; ++i) {
      mpz s; mpz_init(s.get_mpz_t());
      mpz_urandomm(s.get_mpz_t(), state, order); // s in [0, r-1]
      vec.push_back(s);
    }
  };

  std::cout << "Generating scalars...\n";
  gen_scalars(scalars1, N1);
  gen_scalars(scalars2, N2);
  std::cout << "Scalars generated.\n";

  // Benchmark G1 exponentiations
  element_t out1;
  element_init_G1(out1, pairing);
  auto t0 = std::chrono::high_resolution_clock::now();
  for (size_t i = 0; i < N1; ++i) {
    element_pow_mpz(out1, g1.ptr(), scalars1[i].get_mpz_t());
  }
  auto t1 = std::chrono::high_resolution_clock::now();
  double secs_g1 = std::chrono::duration<double>(t1 - t0).count();
  std::cout << "G1: " << N1 << " exponentiations in " << secs_g1 << " s, "
            << (N1 / secs_g1) << " ops/s\n";

  // Benchmark G2 exponentiations
  element_t out2;
  element_init_G2(out2, pairing);
  t0 = std::chrono::high_resolution_clock::now();
  for (size_t i = 0; i < N2; ++i) {
    element_pow_mpz(out2, g2.ptr(), scalars2[i].get_mpz_t());
  }
  t1 = std::chrono::high_resolution_clock::now();
  double secs_g2 = std::chrono::duration<double>(t1 - t0).count();
  std::cout << "G2: " << N2 << " exponentiations in " << secs_g2 << " s, "
            << (N2 / secs_g2) << " ops/s\n";

  // Precompute lists of G1 and G2 elements for pairing to avoid repeated exponentiation inside pairing loop.
//   std::vector<element_t> g1_list;
//   std::vector<element_t> g2_list;
  std::vector<std::unique_ptr<PBCElement>> g1_list;
  std::vector<std::unique_ptr<PBCElement>> g2_list;
  g1_list.reserve(NP);
  g2_list.reserve(NP);
  for (size_t i = 0; i < NP; ++i) {
    auto pe1 = std::make_unique<PBCElement>(pairing, PBCElement::G1);
    element_pow_mpz(pe1->ptr(), g1.ptr(), scalars1[i].get_mpz_t()); // reuse first scalars
    g1_list.push_back(std::move(pe1));
    auto pe2 = std::make_unique<PBCElement>(pairing, PBCElement::G2);
    element_pow_mpz(pe2->ptr(), g2.ptr(), scalars2[i].get_mpz_t()); // reuse first scalars
    g2_list.push_back(std::move(pe2));
  }

  // Benchmark pairings
  auto tp0 = std::chrono::high_resolution_clock::now();
  for (size_t i = 0; i < NP; ++i) {
    element_pairing(tG.ptr(), g1_list[i]->ptr(), g2_list[i]->ptr());
  }
  auto tp1 = std::chrono::high_resolution_clock::now();
  double secs_pair = std::chrono::duration<double>(tp1 - tp0).count();
  std::cout << "GT: " << NP << " pairings in " << secs_pair << " s, "
            << (NP / secs_pair) << " ops/s\n";

  // Cleanup
//   for (auto &s : scalars1) mpz_clear(s.get_mpz_t());
//   for (auto &s : scalars2) mpz_clear(s.get_mpz_t());
  for (auto &e : g1_list) element_clear(e->ptr());
  for (auto &e : g2_list) element_clear(e->ptr());
  element_clear(g1.ptr()); element_clear(g2.ptr()); element_clear(out1); element_clear(out2);
  element_clear(tG.ptr());
  mpz_clear(order);
  mpz_clear(scalar); mpz_clear(tmp);
  gmp_randclear(state);
  pairing_clear(pairing);

  return 0;
}
