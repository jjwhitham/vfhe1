// Performs multi-scalar multiplication (linear combinations) of base points vec
// P = [P_0, ..., P_{N-1}] by set of vecs a = {a_i}_{i=0}^{i=m},
// where a_i = [a_{i,0}, ..., a_{i,N-1}]

// #include <mcl/bn.hpp>
#include "/Users/jw/Projects/mcl/include/mcl/bn.hpp"
#include <iostream>
#include <string>
#include <chrono>
#include <vector>
// TODO include GMP, but may have to wrap in a namespace...

using namespace mcl::bn;

void test_g1_exps(std::vector<std::vector<Fr>>& scalars, size_t N=8192) {
    G1 gen1;
    int gen_seed = 0;
    hashAndMapToG1(gen1, std::string("P_") + std::to_string(gen_seed));
    std::vector<G1> P1;
    // Prepare M scalar vectors (each length N)
    P1.reserve(N);
    auto t0 = std::chrono::high_resolution_clock::now();
    for (size_t i = 0; i < N; ++i) {
        gen1 *= scalars[0][i];
    }
    auto t1 = std::chrono::high_resolution_clock::now();
    double secs = std::chrono::duration<double>(t1 - t0).count();
    std::cout << "G1 exps: " << N << " mults of aG1 in " << secs << " s, "
                        << (N / secs) << " ops/s\n";
}
void test_g2_exps(std::vector<std::vector<Fr>>& scalars, size_t N=8192) {
    G2 gen2;
    int gen_seed = 0;
    hashAndMapToG2(gen2, std::string("P_") + std::to_string(gen_seed));
    std::vector<G2> P2;
    P2.reserve(N);
    // Prepare M scalar vectors (each length N)
    auto t0 = std::chrono::high_resolution_clock::now();
    for (size_t i = 0; i < N; ++i) {
        gen2 *= scalars[0][i];
    }
    auto t1 = std::chrono::high_resolution_clock::now();
    double secs = std::chrono::duration<double>(t1 - t0).count();
    std::cout << "G2 exps: " << N << " mults of aG2 in " << secs << " s, "
                        << (N / secs) << " ops/s\n";
}
void test_gT_exps(std::vector<std::vector<Fr>>& scalars, size_t N=8192) {
    G1 gen1;
    G2 gen2;
    int gen_seed = 0;
    hashAndMapToG1(gen1, std::string("P_") + std::to_string(gen_seed));
    hashAndMapToG2(gen2, std::string("P_") + std::to_string(gen_seed));
    GT gT;
    pairing(gT, gen1, gen2);
    // Prepare M scalar vectors (each length N)
    auto t0 = std::chrono::high_resolution_clock::now();
    for (size_t i = 0; i < N; ++i) {
        GT::pow(gT, gT, scalars[0][i]);
    }
    auto t1 = std::chrono::high_resolution_clock::now();
    double secs = std::chrono::duration<double>(t1 - t0).count();
    std::cout << "GT exps: " << N << " exps of gT^a in " << secs << " s, "
                        << (N / secs) << " ops/s\n";
}
void test_pairing(std::vector<std::vector<Fr>>& scalars, size_t M=300) {
    G1 gen1;
    G2 gen2;
    GT out;
    int gen_seed = 0;
    hashAndMapToG1(gen1, std::string("P_") + std::to_string(gen_seed));
    hashAndMapToG2(gen2, std::string("P_") + std::to_string(gen_seed));
    std::vector<G1> P1;
    std::vector<G2> P2;
    // Prepare M scalar vectors (each length N)
    P1.reserve(M);
    P2.reserve(M);
    for (size_t i = 0; i < M; ++i) {
        P1.push_back(gen1 * scalars[i][0]);
        P2.push_back(gen2 * scalars[i][1]);
    }
    auto t0 = std::chrono::high_resolution_clock::now();
    for (size_t i = 0; i < M; ++i) {
        pairing(out, P1[i], P2[i]);
    }
    auto t1 = std::chrono::high_resolution_clock::now();
    double secs = std::chrono::duration<double>(t1 - t0).count();
    std::cout << "pairing: " << M << " pairings of e(aP, bQ) in " << secs << " s, "
                        << (M / secs) << " ops/s\n";
}
void test_msms(std::vector<std::vector<Fr>>& scalars, std::vector<G1>& T, size_t N=8192, size_t M=300, size_t n_threads=4) {
    // Benchmark optimized multi-scalar (G1::mulVec)
    // mulVec: 140 MSMs of size 8192 in 6.29867 s, 22.2269 ops/s - 16 threads
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
    // mulVecMT: 140 MSMs of size 8192 in 0.878884 s, 159.293 ops/s - 16 threads
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
    // naive: 140 MSMs of size 8192 in 35.3093 s, 3.96496 ops/s  - 16 threads
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
}

std::vector<G1> create_powers_of_t_(G1& gen1, size_t N=8192) {
    // Create T: powers of t multiplied by G1's generator, i.e. (P, tP, (t^2)P, ..., (t^{N-1})P)
    std::vector<G1> T;
    T.reserve(N);
    G1 tpowGen = gen1;
    Fr t = 42;
    for (size_t i = 0; i < N; i++) {
        tpowGen *= t;
        T.push_back(tpowGen);
    }
    return T;
}

std::vector<std::vector<G1>> create_pows_of_pows_of_t_(G1& gen1, size_t N=8192) {
    std::vector<G1> pows_of_t_ = create_powers_of_t_(gen1, N);
    size_t len_pows = 256;
    std::vector<std::vector<G1>> ret;
    ret.reserve(N);
    for (size_t i = 0; i < N; i++) {
        std::vector<G1> pows_of_pow_of_t_;
        pows_of_pow_of_t_.reserve(len_pows);
        pows_of_pow_of_t_.push_back(pows_of_t_.at(i));
        for (size_t j = 1; j < len_pows; j++) {
            pows_of_pow_of_t_.push_back(pows_of_pow_of_t_.at(j - 1) + pows_of_pow_of_t_.at(j - 1));
        }
        ret.push_back(pows_of_pow_of_t_);
    }
    return ret;
}


void test_msms_cache(std::vector<std::vector<Fr>>& scalars, G1& gen1, size_t N=8192, size_t M=30) {
    std::vector<std::vector<G1>> tt = create_pows_of_pows_of_t_(gen1, N);
    std::vector<G1> res;
    res.reserve(M);
    size_t len_pows = tt.size(0).size();
    assert(len_pows == 256);
    for (size_t i = 0; i < scalars.size(); i++) {
        for (size_t j = 0; j < N; j++) {
            G1 msm;
            msm.clear();
            for (size_t k = 0; k < len_pows; k++) {
                // TODO this isn't binary pow yet. Need to read exponent (scalar)
                // in binary and then sum the 1's entries...
                // TODO revise what is meant to happen here...
                msm += scalars.at(i).at(j) * tt.at(j).at(k);
            }
            res.push_back(msm);
        }
    }
}

int main(int argc, char **argv) {
    size_t n_threads = 4;
    std::cout << "n_threads = " << n_threads << "\n";

    // defaults
    size_t N = 8192; // number of base points
    size_t M = 30;  // number of scalar vectors to apply

    if (argc > 1) N = std::stoul(argv[1]);
    if (argc > 2) M = std::stoul(argv[2]);

    // Initialize curve / pairing (sets params for G1/G2/Fr)
    initPairing(BN_SNARK1);
    std::cout << "|Fp| = " << mcl::Fp::getBitSize() << "\n";
    std::cout << "|Fr| = " << mcl::Fr::getBitSize() << "\n";

    // Set a generator
    G1 gen1;
    bool use_random_gen = false;
    if (use_random_gen) {
        Fr s;
        s.setByCSPRNG();
        hashAndMapToG1(gen1, std::string("P_") + s.getStr());
    } else {
        int gen_seed = 0;
        hashAndMapToG1(gen1, std::string("P_") + std::to_string(gen_seed));
    }


    // Prepare M scalar vectors (each length N)
    std::vector<std::vector<Fr>> scalars;
    scalars.reserve(M);
    for (size_t i = 0; i < M; ++i) {
        scalars.emplace_back();
        scalars.back().reserve(N);
        for (size_t j = 0; j < N; ++j) {
            Fr s; s.setByCSPRNG();
            scalars.back().push_back(s);
        }
    }

    std::vector<G1> T = create_powers_of_t_(gen1, N);
    test_g1_exps(scalars, N);
    test_g2_exps(scalars, N);
    test_gT_exps(scalars, N);
    test_pairing(scalars, M);
    test_msms(scalars, T, N, M, n_threads);
    test_msms_cache(scalars, gen1, N, M);
    return 0;
}