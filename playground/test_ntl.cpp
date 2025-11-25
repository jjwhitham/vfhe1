// https://libntl.org/doc/tour-unix.html
//   - see 'After Installing NTL' for compiling/linking/running

#include <NTL/ZZ_pX.h>
#include <NTL/ZZ.h>
#include <assert.h>
#include <iostream>
#include <chrono>
#include <NTL/BasicThreadPool.h>
#include <thread>

using namespace NTL;

void conv_to_nega() {

}

void test_basic_arithmetic() {
    // add
    // sub
    // mul (nega)
}

void test_optimised_conv() {
    // apply all optimisations to negacyclic NTT

    // put x into NTT form
    // for each b in b_vec -> res += NTT(b) * x_NTT
    // iNTT(res)
    const char* q_str = "\
21888242871839275222246405745257275088548364400416034343698204186575808495617";
    ZZ q;
    q = conv<ZZ>(q_str);
    ZZ_p::init(q);

    size_t N = 8192;
    assert(q % (2 * N) == 1);

    // long num_polys = 10000;
    long num_polys = 280;
    Vec<ZZ_pX> polys1;
    Vec<ZZ_pX> polys2;
    polys1.SetLength(num_polys);
    polys2.SetLength(num_polys);
    for (long i = 0; i < num_polys; i++) {
        polys1[i].SetLength(N);
        polys2[i].SetLength(N);
        for (long j = 0; j < (long)N; j++) {
            polys1[i][j] = random_ZZ_p();
            polys2[i][j] = random_ZZ_p();
        }
    }

    long num_threads = AvailableThreads();
    std::cout << "num_threads: " << num_threads << "\n";

    // Vec<ZZ_pX> thread_results;
    // thread_results.SetLength(num_threads);

    // NTT of all poly1s
    long k = 14;
    Vec<FFTRep> polys1_ntt;
    std::cout << "got to here1\n";
    polys1_ntt.SetLength(num_polys);
    std::cout << "got to here2\n";
    for (long i = 0; i < num_polys; i++) {
        std::cout << "got to here3\n";
        FFTRep x(INIT_SIZE, k);
        polys1_ntt[i] = x;
        std::cout << "got to here4\n";
    }

    ZZ_pContext context;
    context.save();
    NTL_EXEC_RANGE(num_polys, first, last)
        context.restore();
        for (long i = first; i < last; i++) {
            ToFFTRep(polys1_ntt[i], polys1[i], k);
        }
    NTL_EXEC_RANGE_END

    // mul + accum of poly1_ntt * NTT(poly2)
    auto start = std::chrono::high_resolution_clock::now();

    NTL_EXEC_RANGE(num_polys, first, last)
        context.restore();
        ZZ_pX accum, tmp;
        accum.SetLength(2 * N);
        tmp.SetLength(2 * N);
        for (long i = first; i < last; i++) {
            FFTRep polys2_ntt(INIT_SIZE, k);
            ToFFTRep(polys2_ntt, polys2[i], k);
            mul(polys1_ntt[i], polys1_ntt[i], polys2_ntt);
            accum += tmp;
        }
    NTL_EXEC_RANGE_END

    auto end = std::chrono::high_resolution_clock::now();
    auto elapsed = std::chrono::duration<double>(end - start).count();
    std::cout << "(mt) " << num_polys << "x poly mults with N=" << N << " took: " << elapsed
              << "(s) (" << num_polys / elapsed << "ops/sec) using " << num_threads << " threads\n";

    ZZ_pX h;
    h.SetLength(2 * N);
    FromFFTRep(h, polys1_ntt[0], 0, 2 * N - 1);
    // // non-destructive versions of the above
    // void NDFromFFTRep(ZZ_pX& x, const FFTRep& y, long lo, long hi, FFTRep& temp);
    // void NDFromFFTRep(ZZ_pX& x, const FFTRep& y, long lo, long hi);
}

void test_optimised_conv1() {
    // apply all optimisations to negacyclic NTT

    // put x into NTT form
    // for each b in b_vec -> res += NTT(b) * x_NTT
    // iNTT(res)
    const char* q_str = "\
21888242871839275222246405745257275088548364400416034343698204186575808495617";
    ZZ q;
    q = conv<ZZ>(q_str);
    ZZ_p::init(q);

    const long n = 4;

    // ZZ_pX poly_mod; // TODO assign const/lead coeffs
    // poly_mod.SetLength(2 * n + 1);
    // // define as per X^N - 1, i.e. cyclic conv (use X^N truncate instead?)
    // poly_mod[0] = -1;
    // poly_mod[2 * n] = 1;
    // ZZ_pXModulus M(poly_mod);
    // ZZ_pXMultiplier G(g, M);

    ZZ_pX f, g;
    f.SetLength(n);
    g.SetLength(n);

    ZZ_pX h;
    h.SetLength(2 * n);

    long k = 3;
    FFTRep f_fft(INIT_SIZE, k);
    ToFFTRep(f_fft, f, k); // k == 3, i.e. 2 * N == 2^k == 2^3

    // TIMING
    FFTRep g_fft(INIT_SIZE, k);
    ToFFTRep(g_fft, g, k); // k == 3, i.e. 2 * N == 2^k == 2^3
    mul(f_fft, f_fft, g_fft);
    // accumulate
    FromFFTRep(h, f_fft, 0, 2 * n - 1);
    // TIMING
}

void test_serialisation() {
    // to and from string
    // to and from mpz_t

    // anything else?
}

void test_features() {
    // set modulus p
    // set some ZZ_p
    // ZZ_p arithmetic

    // create some ZZ_pX
    // add, sub, mul
}
void test_mult_conv() {
    // set x, y, N, p and f = X^N + 1
    // x_conv
    // -------------------------------------------------------------
    // 1. Set the 256-bit prime modulus q
    // -------------------------------------------------------------
    // Example 256-bit prime (you may replace this with your own)
    int x[] = { 8, 7, 10, 5 };
    int y[] = { 2, 6, 9, 13 };
    const char* q_str = "17";
//     const char* q_str = "\
// 21888242871839275222246405745257275088548364400416034343698204186575808495617";
    ZZ q;
    q = conv<ZZ>(q_str);
    ZZ_p::init(q);

    const long n = 4;
    ZZ_pX f, g;
    f.SetLength(n);
    g.SetLength(n);

    for (long i = 0; i < n; i++) {
        // f[i] = random_ZZ_p();
        // g[i] = random_ZZ_p();
        f[i] = x[i];
        g[i] = y[i];
    }
    ZZ_pX poly_mod; // TODO assign const/lead coeffs
    poly_mod.SetLength(2 * n + 1);
    // define as per X^N - 1, i.e. cyclic conv (use X^N truncate instead?)
    poly_mod[0] = -1;
    poly_mod[2 * n] = 1;

    ZZ_pXModulus M(poly_mod);
    ZZ_pXMultiplier G(g, M);
    std::cout << "G = { ";
    for (size_t i = 0; i < (size_t)n; i++)
        std::cout << G.val()[i] << " ";
    std::cout << "}\n";

    ZZ_pX h;
    h.SetLength(2 * n);

    long k = 3;
    FFTRep f_fft(INIT_SIZE, k);
    ToFFTRep(f_fft, f, k); // k == 3, i.e. 2 * N == 2^k == 2^3
    FFTRep g_fft(INIT_SIZE, k);
    ToFFTRep(g_fft, g, k); // k == 3, i.e. 2 * N == 2^k == 2^3
    // FFTRep out_fft(INIT_SIZE, k);
    // ToFFTRep(out_fft, g, k); // k == 3, i.e. 2 * N == 2^k == 2^3
    mul(f_fft, f_fft, g_fft);
    FromFFTRep(h, f_fft, 0, 2 * n - 1);
    // MulMod(h, f, G, M);
    ZZ_pX h1;
    h1.SetLength(n);
    h1 = f * g;
    std::cout << "x * y = { ";
    for (size_t i = 0; i < (size_t)2 * n; i++)
        std::cout << h1[i] << " ";
    std::cout << "}\n";

    // Optional: reduce size
    // h.normalize();

    // print
    // std::cout << "x = { ";
    // for (size_t i = 0; i < (size_t)n; i++)
    //     std::cout << x[i] << " ";
    // std::cout << "}\n";

    // std::cout << "y = { ";
    // for (size_t i = 0; i < (size_t)n; i++)
    //     std::cout << y[i] << " ";
    // std::cout << "}\n";

    std::cout << "xy_conv = { ";
    for (size_t i = 0; i < (size_t)2 * n; i++)
        std::cout << h[i] << " ";
    std::cout << "}\n";

    // degree should be 2(n - 1) == 6
    std::cout << "deg(h) = " << deg(h) << "\n";
}

void test_mult_nega() {
    // set x, y, N, p and f = X^N + 1
    // x_conv
}

void test_ntt_conv() {
    // set x, y, N, p and f = X^N + 1
    // x_conv
    // y_conv
}

void test_ntt_nega() {

}

void test_ntt_peformance() {
    /*
    https://libntl.org/doc/tour-tips.html
    https://libntl.org/doc/tour-time.html
        'multiply degree-1000 int poly with 1000-bit coeffs: 0.00293567'
    */
    // set p to be BN254 group order
    const char* q_str = "\
21888242871839275222246405745257275088548364400416034343698204186575808495617";
    // const char* q_str = "114689";
    ZZ q;
    q = conv<ZZ>(q_str);
    ZZ_p::init(q);
    // assert that q is NTT friendly
    size_t N = 8192;
    assert(q % (2 * N) == 1);
    // create a vec_ZZ_pX of size 280 random polynomials (each poly has 8192 elements)
    size_t num_polys = 280;
    Vec<ZZ_pX> polys;
    polys.SetLength(num_polys);
    for (size_t i = 0; i < num_polys; i++) {
        polys[i].SetLength(N);
        for (size_t j = 0; j < N; j++) {
            polys[i][j] = random_ZZ_p();
        }
    }

    // For 280x iterations, multiply i'th poly with N - 1 - i'th poly
    ZZ_pX result;
    result.SetLength(N);
    auto start = std::chrono::high_resolution_clock::now();
    for (size_t i = 0; i < num_polys; i++) {
        // use num_polys here (was previously using N which caused out-of-range access)
        result = result + polys[i] * polys[num_polys - 1 - i];
    }
    auto end = std::chrono::high_resolution_clock::now();
    auto elapsed = std::chrono::duration<double>(end - start).count();
    std::cout << num_polys << "x poly mults with N=" << N << " took: " << elapsed << "(s) (" << num_polys / elapsed << "ops/sec)\n";

    // TIMING
    // for 1000x iterations, multiply a random pair of polys
    // TIMING
}

void test_ntt_peformance_mt() {
    // https://libntl.org/doc/tour-ex7.html
    const char* q_str = "\
21888242871839275222246405745257275088548364400416034343698204186575808495617";
    ZZ q;
    q = conv<ZZ>(q_str);
    ZZ_p::init(q);

    size_t N = 8192;
    assert(q % (2 * N) == 1);

    long num_polys = 10000;
    Vec<ZZ_pX> polys1;
    Vec<ZZ_pX> polys2;
    polys1.SetLength(num_polys);
    polys2.SetLength(num_polys);
    for (long i = 0; i < num_polys; i++) {
        polys1[i].SetLength(N);
        polys2[i].SetLength(N);
        for (long j = 0; j < (long)N; j++) {
            polys1[i][j] = random_ZZ_p();
            polys2[i][j] = random_ZZ_p();
        }
    }

    long num_threads = AvailableThreads();
    std::cout << "num_threads: " << num_threads << "\n";

    // Vec<ZZ_pX> thread_results;
    // thread_results.SetLength(num_threads);

    auto start = std::chrono::high_resolution_clock::now();

    ZZ_pContext context;
    context.save();
    NTL_EXEC_RANGE(num_polys, first, last)
        fprintf(stderr, "tid %s\n", CurrentThreadID().c_str());
        context.restore();
        ZZ_pX accum, tmp;
        accum.SetLength(N);
        tmp.SetLength(N);
        for (long i = first; i < last; i++) {
            mul(tmp, polys1[i], polys2[i]);
            add(accum, accum, tmp);
            accum += tmp;
        }
            // thread_result += polys[i] * polys[i];
            // thread_result += polys[i] * polys[last - 1 - i];
    NTL_EXEC_RANGE_END

    auto end = std::chrono::high_resolution_clock::now();
    auto elapsed = std::chrono::duration<double>(end - start).count();
    std::cout << "(mt) " << num_polys << "x poly mults with N=" << N << " took: " << elapsed
              << "(s) (" << num_polys / elapsed << "ops/sec) using " << num_threads << " threads\n";
}

void test_ntt_features() {
    // given a poly x and a set of polys Y = {y_0, ..., y_{m-1}}, put x into
    // NTT form and then perform the multiplication with each y_i


    /**************************************************************************\
    from: https://libntl.org/doc/ZZ_pX.cpp.html
                        Modular Arithmetic with Pre-Conditioning

    If you need to do a lot of arithmetic modulo a fixed f, build a
    ZZ_pXModulus F for f.  This pre-computes information about f that
    speeds up subsequent computations.

    It is required that deg(f) > 0 and that LeadCoeff(f) is invertible.

    As an example, the following routine computes the product modulo f of a vector
    of polynomials.

    #include <NTL/ZZ_pX.h>

    void product(ZZ_pX& x, const vec_ZZ_pX& v, const ZZ_pX& f)
    {
    ZZ_pXModulus F(f);
    ZZ_pX res;
    res = 1;
    long i;
    for (i = 0; i < v.length(); i++)
        MulMod(res, res, v[i], F);
    x = res;
    }

    Note that automatic conversions are provided so that a ZZ_pX can
    be used wherever a ZZ_pXModulus is required, and a ZZ_pXModulus
    can be used wherever a ZZ_pX is required.

    \**************************************************************************/


    /**************************************************************************\
    from: https://libntl.org/doc/ZZ_pX.cpp.html
                                    More Pre-Conditioning

    If you need to compute a * b % f for a fixed b, but for many a's, it
    is much more efficient to first build a ZZ_pXMultiplier B for b, and
    then use the MulMod routine below.

    Here is an example that multiplies each element of a vector by a fixed
    polynomial modulo f.

    #include <NTL/ZZ_pX.h>

    void mul(vec_ZZ_pX& v, const ZZ_pX& b, const ZZ_pX& f)
    {
    ZZ_pXModulus F(f);
    ZZ_pXMultiplier B(b, F);
    long i;
    for (i = 0; i < v.length(); i++)
        MulMod(v[i], v[i], B, F);
    }
    \**************************************************************************/
}


void test_parallel() {
}

int main() {
    // test_mult_conv();
    // test_ntt_peformance();
    long nt = 16;
    SetNumThreads(nt);
    // test_ntt_peformance_mt();
    test_optimised_conv();
    return 0;
}