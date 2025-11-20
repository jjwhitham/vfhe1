// test_mcl.cpp
// Benchmark G1/G2 exponentiation and pairing using mcl
// #include "/home/jw/Projects/mcl/include/mcl/bls12_381.hpp"
#include <mcl/bls12_381.hpp>
#include <gmpxx.h>
#include <chrono>
#include <iostream>
#include <vector>
#include <memory>

using namespace mcl::bn;
using mpz = mpz_class;

int main() {
  // initialize pairing for BN curve
  // initPairing() will choose a default BN curve (BN254/BN_P)
  initPairing();

  // generators
  G1 P;
  G2 Q;
  Fr s;

  // set generators randomly (for benchmark purposes)
  P.rand();
  Q.rand();

  // Pre-generate scalars
  const size_t N1 = 1000;
  const size_t N2 = 1000;
  const size_t NP = 100;

  std::vector<Fr> scalars1; scalars1.reserve(N1);
  std::vector<Fr> scalars2; scalars2.reserve(N2);

  gmp_randstate_t state;
  gmp_randinit_default(state);
  gmp_randseed_ui(state, (unsigned long)time(nullptr));

  // get group order as big integer (Fr::getOrder not available portably), use random Fr::setByCSPRNG
  for (size_t i = 0; i < N1; ++i) {
    Fr t; t.setByCSPRNG();
    scalars1.push_back(t);
  }
  for (size_t i = 0; i < N2; ++i) {
    Fr t; t.setByCSPRNG();
    scalars2.push_back(t);
  }

  // Benchmark G1 exponentiations
  G1 out1;
  auto t0 = std::chrono::high_resolution_clock::now();
  for (size_t i = 0; i < N1; ++i) {
    // out1 = scalars1[i] * P
    mul(out1, P, scalars1[i]);
  }
  auto t1 = std::chrono::high_resolution_clock::now();
  double secs_g1 = std::chrono::duration<double>(t1 - t0).count();
  std::cout << "G1: " << N1 << " exponentiations in " << secs_g1 << " s, "
            << (N1 / secs_g1) << " ops/s\n";

  // Benchmark G2 exponentiations
  G2 out2;
  t0 = std::chrono::high_resolution_clock::now();
  for (size_t i = 0; i < N2; ++i) {
    mul(out2, Q, scalars2[i]);
  }
  t1 = std::chrono::high_resolution_clock::now();
  double secs_g2 = std::chrono::duration<double>(t1 - t0).count();
  std::cout << "G2: " << N2 << " exponentiations in " << secs_g2 << " s, "
            << (N2 / secs_g2) << " ops/s\n";

  // Precompute lists of G1 and G2 elements for pairing
  std::vector<std::unique_ptr<G1>> g1_list;
  std::vector<std::unique_ptr<G2>> g2_list;
  g1_list.reserve(NP);
  g2_list.reserve(NP);
  for (size_t i = 0; i < NP; ++i) {
    auto p1 = std::make_unique<G1>();
    auto p2 = std::make_unique<G2>();
    mul(*p1, P, scalars1[i % scalars1.size()]);
    mul(*p2, Q, scalars2[i % scalars2.size()]);
    g1_list.push_back(std::move(p1));
    g2_list.push_back(std::move(p2));
  }

  // Benchmark pairings
  GT e;
  auto tp0 = std::chrono::high_resolution_clock::now();
  for (size_t i = 0; i < NP; ++i) {
    pairing(e, *g1_list[i], *g2_list[i]);
  }
  auto tp1 = std::chrono::high_resolution_clock::now();
  double secs_pair = std::chrono::duration<double>(tp1 - tp0).count();
  std::cout << "GT: " << NP << " pairings in " << secs_pair << " s, "
            << (NP / secs_pair) << " ops/s\n";

  return 0;
}
