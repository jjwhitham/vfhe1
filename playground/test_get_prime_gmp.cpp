// montmul_compare_montgomery.cpp
#include <gmpxx.h>
#include <iostream>
#include <vector>
#include <chrono>
#include <cstdlib>
#include <cstring>

typedef __uint128_t u128;

void find_q(mpz_class& q, mpz_class const& N, mpz_class const& q_init) {
    int i = 0;
    mpz_class k = (q_init - 1) / (2 * N) + 1;
    q = q_init;
    if (mpz_probab_prime_p(q.get_mpz_t(), 25) && (q % (2 * N) == 1))
        return;
    while (i < 100) {
        mpz_nextprime(q.get_mpz_t(), q.get_mpz_t());
        if (q % (2 * N) == 1)
            return;
        k += 1;
        // q = 1 + 2 * N * k; // FIXME drop the 2?
        i += 1;
    }
}

int main(int argc, char* argv[]) {
    if (argc != 3) {
        std::cerr << "Usage: " << argv[0] << " <N> <q_init>\n";
        return 1;
    }
    int N_ = std::atoi(argv[1]);
    if (N_ <= 0) {
        std::cerr << "N must be positive.\n";
        return 1;
    }

    int q_init_ = std::atoi(argv[2]);
    if (q_init_ <= 0) {
        std::cerr << "q_init must be positive.\n";
        return 1;
    }

    mpz_class N = N_;
    mpz_class q_init = q_init_;
    mpz_class q;

    find_q(q, N, q_init);
    // mpz_nextprime(q.get_mpz_t(), q_init.get_mpz_t());
    gmp_printf ("q: %Zd\n", q);
    // std::cout << "q: " << q.get_mpz_t()[0] << "\n";
    // auto end_mont = std::chrono::high_resolution_clock::now();
    // double total_ms_mont = std::chrono::duration<double,std::milli>(end_mont - start_mont).count();


    return 0;
}
