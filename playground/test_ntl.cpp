#include <NTL/ZZ_pX.h>
#include <NTL/ZZ.h>

using namespace NTL;

int main() {
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

    return 0;
}
