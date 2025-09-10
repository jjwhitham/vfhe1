#include "vfhe.h"
#include "shared.h"
#include "enc.h"

void test_rlwe() {
    static constexpr u32 d = N_DECOMP;
    static constexpr u32 power = get_decomp_power();
    static constexpr i128 v = static_cast<i128>(1) << power;
    static constexpr i128 N = N_;
    static constexpr i128 q = FIELD_MODULUS;
    Encryptor enc(v, d, N, q);
    i128 scale = 1000; // 1e3
    i128 message = 42 * scale;
    poly ptx(N);
    ptx.set(0, message); // set the first coefficient to the value
    rlwe ctx = enc.encrypt_rlwe(ptx);
    poly ptx1 = enc.decrypt_rlwe(ctx);
    i128 message_decrypted = ptx1.get(0);
    std::cout << "message decrypted: "
        << i128str(message_decrypted) << "\n";

    /* ### convolution ### */
    rgsw rg = enc.encrypt_rgsw(ptx);
    for (size_t i = 0; i < 2 * d; i++) {
        poly p = enc.decrypt_rlwe(rg.get_rlwe(i));
        std::cout << i128str(p.get_coeff(0)) << " ";
    }
    std::cout << "\n";
    rlwe prod2 = rg * ctx.decompose(v, d, power);

    /* ### negacyclic convolution ### */
    // rg.conv_to_ntt();
    // rlwe prod = rg.convolve(ctx.decompose(v, d, power));
    // std::cout << "yoohoo\n";
    // prod.conv_to_coeff();
    // rlwe prod2 = prod.conv_to_nega(N);

    i128 prod_dec = enc.decrypt_rlwe(prod2).get(1);
    std::cout << "prod decrypted: "
        << i128str(prod_dec) << "\n";
}

void test_decomp() {
    static constexpr u32 d = N_DECOMP;
    static constexpr u32 power = get_decomp_power();
    static constexpr i128 v = static_cast<i128>(1) << power;
    static constexpr i128 N = N_;
    static constexpr i128 q = FIELD_MODULUS;
    Encryptor enc(v, d, N, q);
    i128 scale = 1e12; // 1e12
    // std::cout << "scale: " << i128str(scale) << "\n";
    i128 message = 42 * scale;
    poly ptx(N);
    ptx.set(0, message); // set the first coefficient to the value
    rlwe ctx = enc.encrypt_rlwe(ptx);
    rlwe_decomp ctx_decomp = ctx.decompose(v, d, power);

    // std::cout << "ctx_decomp.size(): " << ctx_decomp.size() << "\n";

    vector_i128 v_pows(d);
    for (size_t i = 0; i < d; i++)
        v_pows.at(i) = pow_(v, i, q);
    // std::cout << "v_pows:\n";
    // print_vector_i128(v_pows);

    poly b(N);
    poly a(N);
    for (size_t i = 0; i < d; i++) {
        b = b + (ctx_decomp.get(i) * v_pows.at(i));
        a = a + (ctx_decomp.get(i + d) * v_pows.at(i));
    }
    rlwe ctx1(N_POLYS_IN_RLWE);
    ctx1.set(0, b);
    ctx1.set(1, a);
    poly ptx1 = enc.decrypt_rlwe(ctx1);
    // for (size_t i = 0; i < 4; i++)
    //     std::cout << "ptx[" << i << "]: " << i128str(ptx1.get(i)) << "\n";
    // std::cout << "message decrypted: " << i128str(ptx1.get(0)) << "\n";
    assert(ptx1.get(0) < 100 - message);
}

int main() {
    test_rlwe();
    test_decomp();
    return 0;
}