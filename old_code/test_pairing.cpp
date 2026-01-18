// #include <pbc/pbc.h>
// #include <assert.h>

// int main() {
//     pbc_param_t param;
//     pairing_t pairing;
//     element_t g1, g2, gt1, gt2, gt3, a, g1a;
//     pbc_param_init_a_gen(param, 160, 512);
//     pairing_init_pbc_param(pairing, param);
//     element_init_G1(g1, pairing);
//     element_init_G2(g2, pairing);
//     element_init_G1(g1a, pairing);
//     element_init_GT(gt1, pairing);
//     element_init_GT(gt2, pairing);
//     element_init_GT(gt3, pairing);
//     element_init_Zr(a, pairing);
//     element_random(g1); element_random(g2); element_random(a);
//     element_pairing(gt1, g1, g2); // gt1 = e(g1, g2)
//     element_pow_zn(g1a, g1, a); // g1a = g1^a
//     element_pow_zn(gt2, gt1, a); // gt2 = gt1^a = e(g1, g2)^a
//     element_pairing(gt3, g1a, g2); // gt3 = e(g1a, g2) = e(g1^a, g2)
//     assert(element_cmp(gt2, gt3) == 0); // assert gt2 == gt3
//     element_clear(g1); element_clear(g2); element_clear(gt1);
//     element_clear(gt2); element_clear(gt3); element_clear(a);
//     element_clear(g1a);
//     pairing_clear(pairing);
//     return 0;
// }

/*
  example.c
  Demonstrates PBC pairing usage with GMP mpz_t:
  - Initialize pairing from parameter string
  - Use mpz_t scalars
  - Get group generators (G1/G2)
  - Point addition, scalar multiplication
  - Compute pairing e : G1 x G2 -> GT
*/

#include <pbc/pbc.h>
#include <gmp.h>
#include <stdio.h>

int main() {
  pairing_t pairing;
  element_t g1, g2;        // generators: G1 and G2
  element_t a1, a2;        // group elements
  element_t t1, t2;        // results
  element_t gt, gt1, gt2;  // target group elements
  mpz_t k;                 // scalar as GMP integer

  // 1) Example parameter string: Type A (symmetric) — replace with BN parameters if available.
  // Type A parameters are of form: "type a\nq ...\nh ...\n"
  // For demo we use built-in generator function pbc_param_init_a_gen to produce params.
  pbc_param_t param;
  // Generate Type A parameters with rbits and qbits (small for demo)
  pbc_param_init_a_gen(param, 256, 256); // rbits=160 (group order size), qbits=512 (field size)
  pairing_init_pbc_param(pairing, param);

  // 2) Initialize elements
  // For Type A pairing G1 == G2 (symmetric); still demonstrate API
  element_init_G1(g1, pairing);
  element_init_G2(g2, pairing); // in symmetric pairings G1==G2 but type differs
  element_init_G1(a1, pairing);
  element_init_G1(a2, pairing);
  element_init_GT(gt, pairing);
  element_init_GT(gt1, pairing);
  element_init_GT(gt2, pairing);

  // 3) Set generators: for Type A we can use random or use hash-to-point
  element_random(g1);
  element_set(g2, g1); // symmetric example: set g2 = g1

  // 4) Demonstrate addition
  element_set(a1, g1);
  element_add(a2, a1, g1); // a2 = a1 + g1  (i.e., 2*g1)
  element_printf("a2 = a1 + g1: %B\n", a2);

  // 5) Scalar multiplication using GMP mpz_t
  mpz_init_set_str(k, "12345678901234567890", 10); // sample scalar
  // PBC supports element_pow_mpz to exponentiate in GT or scalar multiply in G1 (if appropriate)
  element_pow_mpz(a1, g1, k); // a1 = g1^k  (in additive notation for G1 this is scalar * g1)
  element_printf("a1 = k * g1: %B\n", a1);

  // 6) Pairing computation
  element_pairing(gt, a1, g2); // gt = e(a1, g2)
  element_printf("pairing e(a1,g2): %B\n", gt);

  // 7) Check bilinearity: e(g1^k, g2) == e(g1, g2)^k
  element_pairing(gt1, a1, g2); // e(g1^k, g2)
  element_pairing(gt2, g1, g2); // e(g1, g2)
  element_pow_mpz(gt2, gt2, k); // gt2 = e(g1,g2)^k
  int equal = element_cmp(gt1, gt2) == 0;
  printf("Bilinearity holds? %s\n", equal ? "yes" : "no");

  // 8) Cleanup
  mpz_clear(k);
  element_clear(g1);
  element_clear(g2);
  element_clear(a1);
  element_clear(a2);
  element_clear(gt);
  element_clear(gt1);
  element_clear(gt2);
  pairing_clear(pairing);
  pbc_param_clear(param);

  return 0;
}
