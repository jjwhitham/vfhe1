// https://libntl.org/doc/tour-unix.html
//   - see 'After Installing NTL' for compiling/linking/running

#include <NTL/ZZ_pX.h>
#include <NTL/ZZ.h>

using namespace NTL;

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
    ZZ_pX h;
    h.SetLength(2 * n);

    long k = 3;
    FFTRep f_fft(INIT_SIZE, 2 * n);
    ToFFTRep(f_fft, f, k); // k == 3, i.e. 2 * N == 2^k == 2^3
    FFTRep g_fft(INIT_SIZE, 2 * n);
    ToFFTRep(g_fft, g, k); // k == 3, i.e. 2 * N == 2^k == 2^3
    // FFTRep out_fft(INIT_SIZE, 2 * n);
    // ToFFTRep(out_fft, g, k); // k == 3, i.e. 2 * N == 2^k == 2^3
    mul(f_fft, f_fft, g_fft);
    FromFFTRep(h, f_fft, 0, 2 * n - 1);
    // MulMod(h, f, G, M);
    // ZZ_pX h = f * g;

    // Optional: reduce size
    // h.normalize();

    // degree should be 2(n - 1) == 6
    std::cout << "deg(h) = " << deg(h) << "\n";

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

    std::cout << "G = { ";
    for (size_t i = 0; i < (size_t)n; i++)
        std::cout << G.val()[i] << " ";
    std::cout << "}\n";
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
    // assert that p is NTT friendly
    // create 1000x random polynomials with Zp coeff vec of length 4096

    // TIMING
    // for 1000x iterations, perform NTT on poly
    // TIMING


    // TIMING
    // for 1000x iterations, multiply a random pair of polys
    // TIMING
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
    // https://libntl.org/doc/tour-ex7.html
}

void test() {
    // -------------------------------------------------------------
    // 1. Set the 256-bit prime modulus q
    // -------------------------------------------------------------
    // Example 256-bit prime (you may replace this with your own)
    const char* q_str =
        "115792089237316195423570985008687907853269984665640564039457"
        "584007908834671";   // = 2^256 - 189 (just an example)

    ZZ q;
    q = conv<ZZ>(q_str);

    ZZ_p::init(q);   // GF(q)

    // -------------------------------------------------------------
    // 2. Construct two polynomials f(x) and g(x)
    //    Degree <= 4096
    // -------------------------------------------------------------
    const long n = 4096;

    ZZ_pX f, g;

    f.SetLength(n + 1);
    g.SetLength(n + 1);

    // Fill with random coefficients mod q
    for (long i = 0; i <= n; i++) {
        f[i] = random_ZZ_p();
        g[i] = random_ZZ_p();
    }

    // -------------------------------------------------------------
    // 3. Multiply h = f * g   (mod q)
    // -------------------------------------------------------------
    ZZ_pX h = f * g;

    // Optional: reduce size
    h.normalize();

    // -------------------------------------------------------------
    // 4. Output degree of result
    // -------------------------------------------------------------
    std::cout << "deg(h) = " << deg(h) << "\n";  // should be 8192
}

int main() {
    test_mult_conv();
    return 0;
}