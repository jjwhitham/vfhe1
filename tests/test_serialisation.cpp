#include "serialisation.hpp"
#include <gmpxx.h>
#include <NTL/ZZ.h>
#include <NTL/ZZ_p.h>
#include <gmpxx.h>
#include "/home/jw/Projects/mcl/include/mcl/bn.hpp"
// #include <mcl/bn.hpp>
#include <boost/exception/diagnostic_information.hpp>

void ntl_to_mcl_with_binary_serialisation() {
    // inline void BytesFromZZ(unsigned char *p, const ZZ& a, long n)
// 30644e72e131a029b85045b68181585d2833e84879b9709143e1f593f0000001 (q)

/*
00|11_0000 0110_0100 0100_1110 0111_0010 1110_0001 0011_0001 1010_0000 0010_1001
   10111000 01010000 01000101 10110110 10000001 10000001 01011000 0101_1101
   00101000 00110011 11101000 01001000 01111001 10111001 01110000 10010001
   01000011 11100001 11110101 10010011 11110000 00000000 00000000 00000000 (q-1)
*/

// 1100000110010001001110011100101110000100110001101000000010100110111000010100000100010110110110100000011000000101011000010111010010100000110011111010000100100001111001101110010111000010010001
// 01000011 11100001 11110101 10010011 11110000 00000000 00000000 00000000 (q-1)
    const char* q_chars = "\
21888242871839275222246405745257275088548364400416034343698204186575808495617"; // q

    // Serialise NTL to a byte buffer
    ZZ q;
    q = conv<ZZ>(q_chars);
    ZZ qneg1 = q - 1; // XXX meant to error?
    long n_bytes = NumBytes(qneg1);
    assert(n_bytes == 32); // 256 bits (2^8) == 2^5 Bytes
    BytesFromZZ(BUF, qneg1, n_bytes);
    bool print_reversed = true; // bytes are stored least significant first
    print_bytes1(BUF, n_bytes, print_reversed);
    // print_bytes_detailed(BUF, n_bytes);

    // De-serialise MCL from the buffer
    mcl::Fr y;
    buff_to_mcl(y, n_bytes, BUF);
    memset(BUF, 0, N_BYTES_256_BITS);

    // Create another MCL from the string and assert equality
    std::ostringstream qneg1_oss;
    qneg1_oss << qneg1;
    std::string qneg1_str = qneg1_oss.str();
    mcl::Fr x;
    x.setStr(qneg1_str, 10);
    // std::cout << x.getStr(10) << " - x.getStr(10)\n";
    assert(x == y);

    // Serialise MCL and print bytes
    uint8_t buf1[N_BYTES_256_BITS] = { 0 };
    size_t n_bytes_return;
    mcl_to_buff(buf1, &n_bytes_return, x);
    assert(n_bytes_return == (size_t)n_bytes);
    print_bytes1(buf1, n_bytes, print_reversed);
    // print_bytes_detailed(buf1, n_bytes);

    ZZ z = ZZFromBytes(buf1, n_bytes);
    memset(buf1, 0, N_BYTES_256_BITS);
    assert(z == qneg1);
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
    mcl::Fr Fr_cpp_style = ZZ_p_to_new_Fr(qneg1);

    std::cout << qneg1_str << "\n";
    std::cout << Fr_c_style.getStr(10) << "\n";
    std::cout << Fr_cpp_style.getStr(10) << "\n";

    assert(qneg1_str == Fr_c_style.getStr(10));
    assert(qneg1_str == Fr_cpp_style.getStr(10));
}
/* Set rop from an array of word data at op.
The parameters specify the format of the data. count many words are read, each size bytes.
order can be 1 for most significant word first or -1 for least significant first. Within each
word endian can be 1 for most significant byte first, -1 for least significant first, or 0 for the
native endianness of the host CPU. The most significant nails bits of each word are skipped,
this can be 0 to use the full words.
There is no sign taken from the data, rop will simply be a positive integer. An application
can handle any sign itself, and apply it for instance with mpz_neg.
There are no data alignment restrictions on op, any address is allowed.
Here’s an example converting an array of unsigned long data, most significant element first,
and host byte order within each value. */
// void mpz_import1(
//     mpz_t rop, // output (return operand)
//     size_t count, // n_words (32)
//     int order, // word order (-1 least significant first)
//     size_t size, // n_bytes_per_word (1)
//     int endian, // little endian (-1), or host (0)
//     size_t nails, // ?? kind of padding?
//     const void *op) {} // operand

// // `the return value is the destination used, either rop or the allocated block.'
// void* mpz_export1(
//     void *rop, // return operand
//     size_t *countp, // n_words written (passing NULL is ok too)
//     int order, // -1 for least significant word first
//     size_t size, // n_bytes_per_word (1 - 1 * sizeof uint8_t)
//     int endian, // -1 for little endian, 0 for host
//     size_t nails, // padding of the most significant end of each word
//     const mpz_t op) {} // operand

void test_gmp_to_and_from_binary() {
    // TODO mpz from string
    // const mpz_class from_string("21888242871839275222246405745257275088548364400416034343698204186575808495617");
    std::string q_str("18014398509506561");
    const mpz_class from_string(q_str);
    size_t count_out = 42; // a deterministic value we're not to be returned
    // TODO mpz export
    mpz_export(BUF, &count_out, ORDER, WORD_SIZE, ENDIAN, NAILS, from_string.get_mpz_t());
    // assert(count_out == 32);
    std::cout << "Wrote " << count_out << \
        " Bytes to buf for q = " + q_str + " (55 bit prime)\n";
    bool print_reversed = (ENDIAN == 1) ? false : true;
    print_bytes1(BUF, count_out, print_reversed);

    // TODO mpz import
    // init from_bytes to something 256 bits long
    mpz_class from_bytes = 4 * from_string; // 256 bits (from_string is 254)
    mpz_import(from_bytes.get_mpz_t(), N_BYTES_256_BITS, ORDER, WORD_SIZE, ENDIAN, NAILS, BUF);

    // TODO assert equality
    assert(from_bytes == from_string);

    // clear BUF
    memset((void*)BUF, 0, N_BYTES_256_BITS);
}


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


int main() {
    // initPairing(mcl::BN160);     // 2^7 adicity
    // initPairing(mcl::BN_P256);   // 2^2
    initPairing(mcl::BN_SNARK1); // 2^28
    // initPairing(mcl::BLS12_381); // 2^32

    std::string mod = mcl::Fr::getModulo();

    // std::cout << "\nmodulus = " << mod << " (BN160)\n\n";
    // std::cout << "\nmodulus = " << mod << " (BN_P256)\n";
    std::cout << "\nmodulus = " << mod << " (BN_SNARK1)\n\n";
    // std::cout << "\nmodulus = " << mod << " (BLS12_381)\n\n";
    // ntl_to_mcl_with_string();
    // ntl_to_mcl_with_mpz();
    // ntl_to_mcl_with_binary_serialisation();
    // test_ZZ_p_to_Fr();
    test_gmp_to_and_from_binary();
    return 0;
}

