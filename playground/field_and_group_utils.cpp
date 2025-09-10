// #include <gmpxx.h>

// typedef __uint128_t u128;

// void find_q(mpz_class& q, mpz_class const& N, mpz_class const& q_init) {
//     int i = 0;
//     mpz_class k = (q_init - 1) / (2 * N) + 1;
//     q = q_init;
//     if (mpz_probab_prime_p(q.get_mpz_t(), 25) && (q % (2 * N) == 1))
//         return;
//     while (i < 100) {
//         mpz_nextprime(q.get_mpz_t(), q.get_mpz_t());
//         if (q % (2 * N) == 1)
//             return;
//         k += 1;
//         // q = 1 + 2 * N * k; // FIXME drop the 2?
//         i += 1;
//     }
// }

// // Find prime q: q = 1 + 2*N*k, k >= 1, q % (2*N) == 1, q is prime
// i128 find_q(i128 N, i128 q_init) {
//     int i = 0;
//     i128 k = ((q_init - 1) / (2 * N)) + 1;
//     i128 q = 1 + 2 * N * k;
//     while (i < 100) {
//         if ((q % (2 * N)) == 1 && is_prime(q)) {
//             return q;
//         }
//         i++;
//         k++;
//         q = 1 + 2 * N * k;
//     }
//     // If not found, return -1 or throw
//     throw std::runtime_error("No suitable prime q found.");
// }

// // Generation of a cyclic group of prime order
// // Find prime p: p = k * q + 1, k >= 2
// // Given a prime q, find another prime p such that q | (p-1)
// std::pair<i128, i128> get_parent_group(i128 q) {
//     if (!is_prime(q)) {
//         throw std::invalid_argument("q must be a prime number.");
//     }
//     i128 k = 2;
//     while (!is_prime(k * q + 1)) {
//         k++;
//         if (k == 1000000) {
//             throw std::runtime_error("No suitable prime found for the given order.");
//         }
//     }
//     i128 p = k * q + 1;
//     return {p, k};
// }

// // The cyclic group Z/pZ* (of order p - 1) contains a cyclic subgroup of order q.
// // Find a generator of this cyclic subgroup.
// i128 get_generator_of_prime_order_group(i128 p, i128 q, i128 k) {
//     // Randomly pick a number from [2, q-1] and check if it is a generator of the cyclic subgroup of order q.
//     for (i128 h = 2; h < q; ++h) {
//         i128 g = mod_pow(h, k, p);
//         if (g != 1) {
//             // Candidate is a generator of order q
//             if (mod_pow(g, q, p) == 1) {
//                 return g;
//             }
//         }
//     }
//     throw std::runtime_error("No generator found for the given prime order.");
// }

// // For a given finite field of size q (prime), generate a cyclic group of prime order q.
// // q is a prime number ~= 2**54
// // Example: q = (5 * 7 * 62829235873) * 2^13 + 1 (form: m * 2^k + 1 for NTT)
// // q = 18014398509506561
// // q = 17
// // q = 1073750017
// std::tuple<i128, i128, i128> generate_field_and_group_params() {
//     i128 q = find_q(1 << 12, 1 << 12);
//     auto [p, k] = get_parent_group(q);
//     i128 g = get_generator_of_prime_order_group(p, q, k);
//     std::cout << "p: " << p << std::endl;
//     std::cout << "q: " << q << std::endl;
//     std::cout << "k: " << k << std::endl;
//     std::cout << "g: " << g << std::endl << std::endl;
//     // std::cout << "Generated field with prime " << q << " and cyclic group of order " << p-1 << " with generator " << g << "." << std::endl;
//     return {g, q, p};
// }