#include <NTL/ZZ.h>
#include <NTL/ZZ_p.h>
// #include <gmpxx.h>
#include <mcl/bn.hpp>
#include <boost/exception/diagnostic_information.hpp>

using namespace NTL;
// using namespace mcl;

void ntl_to_mcl_with_string() {
    const char* q_chars = "\
21888242871839275222246405745257275088548364400416034343698204186575808495617";
    ZZ q;
    q = conv<ZZ>(q_chars);
    ZZ_p::init(q);
    ZZ_p qneg1 = conv<ZZ_p>(q - 1);
    // ZZ z = rep(q_neg_1); // extract the underlying integer

    std::ostringstream qneg1_oss;
    qneg1_oss << qneg1;
    std::string qneg1_str = qneg1_oss.str();


    mcl::Fr y;
    // XXX this version of setStr throws
    // C++ setStr: void setStr(const std::string& str, int ioMode = 0)
    // y.setStr(qneg1_str);
    std::cout << "\n\nC++ setStr: Setting y to q - 1 ...\n";
    y.setStr(qneg1_str, 10);
    std::cout << y.getStr() << " - y.getStr()\n\n";

    // Try and set y to the size of the modulus
    std::cout << "C++ setStr: Attempting to set y to q ...\n";
    try {
        std::string q_str = std::string(q_chars);
        std::cout << q_str << " - q_str\n";
        y.setStr(q_str, 10);
    } catch(...) {
        // https://stackoverflow.com/a/24142104
        // std::exception_ptr p = std::current_exception();
        // std::clog <<(p ? p.__cxa_exception_type()->name() : "null") << std::endl;
        std::cout << "\nCaught an exception:\n";
        std::clog << boost::current_exception_diagnostic_information() << std::endl;
    }
    std::cout << "After caught exception:\n";
    std::cout << y.getStr() << " - y.getStr()\n\n\n";

    std::cout << "C setStr: Setting y to q - 1 ...\n";
    // XXX this might be faster, as it doesn't throw?
    // C setStr: void setStr(bool *pb, const char *str, int ioMode = 0)
    bool pb;
    // y.setStr(&pb, qneg1_str.c_str());
    y.setStr(&pb, qneg1_str.c_str(), 10);
    assert(pb == true);
    std::cout << y.getStr() << " - y.getStr()\n";
    std::cout << "pb: " << pb << "\n\n";

    std::cout << "C setStr: Attempting to set y to q ...\n";
    std::string mod = mcl::Fr::getModulo();
    std::string q_str = std::string(q_chars);
    assert(mod == q_str);
    y.setStr(&pb, mod.c_str(), 10);
    std::cout << "After failed attempt:\n";
    std::cout << y.getStr() << " - y.getStr()\n";
    std::cout << "pb: " << pb << "\n";
}

void ntl_to_mcl_with_mpz() {
}

void print_bytes(uint8_t *buf, size_t len) {
    for (size_t i = 0; i < len; i++) {
        printf("%02x ", buf[i]);
    }
    printf("\n");
}
void print_byte_in_binary(uint8_t _char, char sep) {
    int n_bits_byte = 8;
    for (int i = n_bits_byte - 1; i >= 0; i--) {
        // printf("%d", (first_char >> (n_bits_byte - 1 - i)) & 1);
        printf("%d", (_char >> i) & 1);
        if (i == 4) {
            printf("_");
        }
    }
    printf("%c", sep);
}

void print_byte_in_hex(const uint8_t *_char) {
    uint8_t upper_nibble = (uint8_t)*_char >> 4;
    uint8_t lower_nibble = *_char & 0b00001111;
    printf("0x%x%x", upper_nibble, lower_nibble);
}

void print_bytes_detailed(const uint8_t *_char, const int n) {
    for (int i = 0; i < n; i++) {
        printf("add: %p, val: 0b", _char);
        char sep = '\0';
        print_byte_in_binary(*_char, sep);
        printf(" | ");
        print_byte_in_hex(_char);
        // if (i % 4 == 3)
        printf(", char: %c", *_char);
        printf("\n");
        _char++;
    }
    printf("\n");
}
void ntl_to_mcl_with_binary_serialisation() {
    // inline void BytesFromZZ(unsigned char *p, const ZZ& a, long n)

/*
00|11_0000 0110_0100 0100_1110 0111_0010 1110_0001 0011_0001 1010_0000 0010_1001
   10111000 01010000 01000101 10110110 10000001 10000001 01011000 0101_1101
   00101000 00110011 11101000 01001000 01111001 10111001 01110000 10010001
   01000011 11100001 11110101 10010011 11110000 00000000 00000000 00000000
*/

// 1100000110010001001110011100101110000100110001101000000010100110111000010100000100010110110110100000011000000101011000010111010010100000110011111010000100100001111001101110010111000010010001
// 01000011 11100001 11110101 10010011 11110000 00000000 00000000 00000000
    const char* q_chars = "\
21888242871839275222246405745257275088548364400416034343698204186575808495617";
    ZZ q;
    q = conv<ZZ>(q_chars);

    ZZ qneg1 = q - 1; // XXX meant to error?
    long n_bytes = NumBytes(qneg1);
    assert(n_bytes == 32); // 256 bits (2^8) == 2^5 Bytes
    unsigned char *buf = (unsigned char *)malloc(sizeof *buf * n_bytes);
    BytesFromZZ(buf, qneg1, n_bytes);
    print_bytes(buf, n_bytes);
    // print_bytes_detailed(buf, n_bytes);

    std::ostringstream qneg1_oss;
    qneg1_oss << qneg1;
    std::string qneg1_str = qneg1_oss.str();
    mcl::Fr x;
    x.setStr(qneg1_str, 10);
    std::cout << x.getStr(10) << " - x.getStr(10)\n";

    uint8_t *buf1 = (uint8_t *)malloc(sizeof *buf1 * n_bytes);
    size_t n_bytes_return = x.getLittleEndian(buf1, n_bytes);
    assert(n_bytes_return == (size_t)n_bytes);
    print_bytes(buf1, n_bytes);
    // print_bytes_detailed(buf1, n_bytes);
    mcl::Fr y;
    y.setLittleEndianMod(buf, n_bytes);
    assert(x == y);
}

// a version that's better if all Fr's are allocated ahead of time
void ZZ_p_to_Fr(mcl::Fr& out, const ZZ_p& in) {
    // temp buffer on stack
    const size_t N_BYTES_256 = 32;
    uint8_t buf[N_BYTES_256] = { 0 }; // TODO zero terminate? zero init?
    size_t n_bytes_in = NumBytes(conv<ZZ>(in));
    assert(n_bytes_in == N_BYTES_256); // XXX this should fail
    BytesFromZZ(buf, conv<ZZ>(in), n_bytes_in); // TODO temporary var for conv
    out.setLittleEndianMod(buf, n_bytes_in);
}
// given a ZZ_p, serialise to a buffer on the stack and return an mcl::Fr
mcl::Fr ZZ_p_to_Fr(const ZZ_p& in) {
    mcl::Fr out;
    ZZ_p_to_Fr(out, in);
    return out;
}

void test_ZZ_p_to_Fr() {
    const char* q_chars = "\
21888242871839275222246405745257275088548364400416034343698204186575808495617";
    ZZ q;
    q = conv<ZZ>(q_chars);
    ZZ_p::init(q);
    ZZ_p qneg1 = conv<ZZ_p>(q - 1);

    std::ostringstream oss;
    oss << qneg1;
    std::string qneg1_str = oss.str();

    mcl::Fr Fr_c_style;
    ZZ_p_to_Fr(Fr_c_style, qneg1);
    mcl::Fr Fr_cpp_style = ZZ_p_to_Fr(qneg1);

    std::cout << qneg1_str << "\n";
    std::cout << Fr_c_style.getStr(10) << "\n";
    std::cout << Fr_cpp_style.getStr(10) << "\n";

    assert(qneg1_str == Fr_c_style.getStr(10));
    assert(qneg1_str == Fr_cpp_style.getStr(10));
}

int main() {
    // initPairing(mcl::BN160);     // 2^7 adicity
    // initPairing(mcl::BN_P256);   // 2^2
    // initPairing(mcl::BN_SNARK1); // 2^28
    initPairing(mcl::BLS12_381); // 2^32

    std::string mod = mcl::Fr::getModulo();

    // std::cout << "\nmodulus = " << mod << " (BN160)\n\n";
    // std::cout << "\nmodulus = " << mod << " (BN_P256)\n";
    // std::cout << "\nmodulus = " << mod << " (BN_SNARK1)\n\n";
    std::cout << "\nmodulus = " << mod << " (BLS12_381)\n\n";
    // ntl_to_mcl_with_string();
    // ntl_to_mcl_with_mpz();
    // ntl_to_mcl_with_binary_serialisation();
    test_ZZ_p_to_Fr();
    return 0;
}