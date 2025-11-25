// Performs multi-scalar multiplication (linear combinations) of base points vec
// P = [P_0, ..., P_{N-1}] by set of vecs a = {a_i}_{i=0}^{i=m},
// where a_i = [a_{i,0}, ..., a_{i,N-1}]

// #include <gmpxx.h>
// typedef mpz_class mpz;
// using mpz = mpz_class;
#include <mcl/bn.hpp>
#include <iostream>
#include <string>
#include <chrono>
#include <vector>
// TODO include GMP, but may have to wrap in a namespace...

using namespace mcl::bn;

int main(int argc, char **argv) {
	size_t n_threads = 16;
    std::cout << "n_threads = " << n_threads << "\n";

	// defaults
	size_t N = 8192; // number of base points
	size_t M = 40;  // number of scalar vectors to apply

	if (argc > 1) N = std::stoul(argv[1]);
	if (argc > 2) M = std::stoul(argv[2]);

	// Initialize curve / pairing (sets params for G1/G2/Fr)
	initPairing();
    std::cout << "|Fp| = " << mcl::Fp::getBitSize() << "\n";
    std::cout << "|Fr| = " << mcl::Fr::getBitSize() << "\n";
	// Prepare base points P
	std::vector<G1> P; P.reserve(N);
	for (size_t i = 0; i < N; ++i) {
		G1 p;
		hashAndMapToG1(p, std::string("P_") + std::to_string(i));
		P.push_back(p);
	}

	// Create T: powers of t multiplied by G1's generator, i.e. (tP, (t^2)P, ..., (t^{N-1})P)
	Fr t = 42;
	Fr t_pow = 1;
	std::vector<G1> T;
	T.reserve(N);
	G1 Generator;
	bool use_random_gen = false;
	if (use_random_gen) {
		Fr s;
		s.setByCSPRNG();
		hashAndMapToG1(Generator, std::string("P_") + s.getStr());
	} else {
		int gen_seed = 0;
		hashAndMapToG1(Generator, std::string("P_") + std::to_string(gen_seed));
	}
	G1 tpowGen = Generator;
	for (size_t i = 0; i < N; i++) {
		// procedural style
		G1::mul(tpowGen, tpowGen, t);
		// alternative declarative/functional style
		tpowGen *= t;

		T.push_back(tpowGen);
	}


	// Prepare M scalar vectors (each length N)
	std::vector<std::vector<Fr>> scalars; scalars.reserve(M);
	for (size_t i = 0; i < M; ++i) {
		scalars.emplace_back();
		scalars.back().reserve(N);
		for (size_t j = 0; j < N; ++j) {
			Fr s; s.setByCSPRNG();
			std::cout << s.getStr(16) << "\n";
			scalars.back().push_back(s);
		}
	}

	// Benchmark optimized multi-scalar (G1::mulVec)
	{
		G1 out;
		auto t0 = std::chrono::high_resolution_clock::now();
		for (size_t i = 0; i < M; ++i) {
			G1::mulVec(out, T.data(), scalars[i].data(), N);
		}
		auto t1 = std::chrono::high_resolution_clock::now();
		double secs = std::chrono::duration<double>(t1 - t0).count();
		std::cout << "mulVec: " << M << " MSMs of size " << N << " in " << secs << " s, "
							<< (M / secs) << " ops/s\n";
	}

	// Benchmark multi-threaded mulVecMT (auto threads)
	{
		G1 out;
		auto t0 = std::chrono::high_resolution_clock::now();
		for (size_t i = 0; i < M; ++i) {
			G1::mulVecMT(out, T.data(), scalars[i].data(), N, n_threads);
		}
		auto t1 = std::chrono::high_resolution_clock::now();
		double secs = std::chrono::duration<double>(t1 - t0).count();
		std::cout << "mulVecMT: " << M << " MSMs of size " << N << " in " << secs << " s, "
							<< (M / secs) << " ops/s\n";
	}

	// Benchmark naive approach: per-scalar multiplication then add
	{
		G1 out;
		auto t0 = std::chrono::high_resolution_clock::now();
		for (size_t i = 0; i < M; ++i) {
			out.clear();
			for (size_t j = 0; j < N; ++j) {
				G1 tmp;
				G1::mul(tmp, T[j], scalars[i][j]);
				out += tmp;
			}
		}
		auto t1 = std::chrono::high_resolution_clock::now();
		double secs = std::chrono::duration<double>(t1 - t0).count();
		std::cout << "naive: " << M << " MSMs of size " << N << " in " << secs << " s, "
							<< (M / secs) << " ops/s\n";
	}

	return 0;
}