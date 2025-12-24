// TODO update size_t* n_read etc. to size_t&
// TODO should the buff_to_xyz funcs do the memset, or should the caller?
#pragma once

// #include "shared.h"
#include <NTL/ZZ.h>
#include <NTL/ZZ_p.h>
#include <gmpxx.h>
#include "/home/jw/Projects/mcl/include/mcl/bn.hpp"
// #include <mcl/bn.hpp>
#include <boost/exception/diagnostic_information.hpp>

static constexpr size_t _256_BITS = 256;
static constexpr size_t N_BYTES_256_BITS = _256_BITS / 8;
static uint8_t BUF[N_BYTES_256_BITS] = { 0 };
static constexpr int ORDER = -1; // least significant word ordering
static constexpr int ENDIAN = -1; // little endian
static constexpr size_t WORD_SIZE = 1 * sizeof (uint8_t); // one Byte words
static constexpr size_t NAILS = 0; // no padding

using namespace NTL;
// using namespace mcl;

void print_bytes1(uint8_t *buf, size_t len, bool reverse=false) {
    if (reverse) {
        for (size_t i = len; i >= 1; i--)
            printf("%02x ", buf[i - 1]);
    } else {
        for (size_t i = 0; i < len; i++)
            printf("%02x ", buf[i]);
    }
    printf("\n");
}
void print_byte_in_binary1(uint8_t _char, char sep) {
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

void print_byte_in_hex1(const uint8_t *_char) {
    uint8_t upper_nibble = (uint8_t)*_char >> 4;
    uint8_t lower_nibble = *_char & 0b00001111;
    printf("0x%x%x", upper_nibble, lower_nibble);
}

void print_bytes_detailed1(const uint8_t *_char, const int n) {
    for (int i = 0; i < n; i++) {
        printf("add: %p, val: 0b", _char);
        char sep = '\0';
        print_byte_in_binary1(*_char, sep);
        printf(" | ");
        print_byte_in_hex1(_char);
        // if (i % 4 == 3)
        printf(", char: %c", *_char);
        printf("\n");
        _char++;
    }
    printf("\n");
}

char* print_to_string_mpz1(const mpz_class& m) {
    int size = 10000;
    char* buf = new char[size];
    buf[size - 1] = '\0';
    int ret = gmp_sprintf(buf, "%Zd", m.get_mpz_t());
    if (ret > size - 1 || ret < 0) {
        throw std::runtime_error("Buffer overflow in print_to_string_mpz");
    }
    return buf;
}


/******************************************************************************/
/**************************** SERIALISATION ***********************************/
/******************************************************************************/
// XXX Only abs value is written. Use mpz_sgn to store sign if required
void mpz_to_buff(uint8_t* buf, size_t* n_write, const mpz_class& in) {
    // mpz_sizeinbase(in.get_mpz_t(), base); // exact for base-2, else exact[+1]
    // mpz_size(in.get_mpz_t()); // n_limbs
    uint8_t* ret_ptr = nullptr;
    ret_ptr = (uint8_t*)mpz_export((void*)buf, n_write, ORDER, WORD_SIZE, ENDIAN, NAILS, in.get_mpz_t());

    // if (ret_ptr == nullptr)
    //     std::cout << "mpz_to_buff: ret_ptr = nullptr\n";
    // else if (ret_ptr == BUF)
    //     std::cout << "mpz_to_buff: ret_ptr points to BUF\n";
    // else {
    //     std::cout << "mpz_to_buff: ret_ptr points to something else (GMP allocator?), ret_ptr = " << ret_ptr << "\n";
    //     throw std::runtime_error("Aborting...\n");
    // }
    // std::cout << *n_write << " Bytes written to the buffer\n";
}

void ntl_to_buff(uint8_t* buf, long* n_write, const ZZ& in) {
    *n_write = NumBytes(in);
    BytesFromZZ(buf, in, *n_write);
}

void mcl_to_buff(uint8_t* buf, size_t* n_write, const mcl::Fr& in) {
    // Ask it to write the full 32 Bytes and see what it comes back with
    *n_write = in.getLittleEndian(buf, N_BYTES_256_BITS);
    // 0 means error
    if (*n_write == 0)
        std::cout << "\nERROR: mcl_to_buff: getLittleEndian returned 0...\n\n";
    // if in == 0, then buf[0] = 0 and return = 1
    if (in == 0) {
        assert(buf[0] == '\0');
        assert(*n_write == 1);
    }
    // otherwise return should be the number of Bytes written
    if (*n_write != N_BYTES_256_BITS)
        std::cout << "\nINFO: mcl_to_buff: getLittleEndian returned " << *n_write << "...\n\n";
    assert(*n_write <= N_BYTES_256_BITS);
    assert(*n_write == N_BYTES_256_BITS);
}

/******************************************************************************/
/**************************** DE-SERIALISATION ********************************/
/******************************************************************************/
// XXX Only absolute value store on buffer. If sign is required, then mpz_sgn
// to access the sign, store it, then use mpz_neg to apply sign on read
void buff_to_mpz(mpz_class& out, size_t n_read, const uint8_t* buf) {
    mpz_import(out.get_mpz_t(), n_read, ORDER, WORD_SIZE, ENDIAN, NAILS, buf);
}

void buff_to_ntl(ZZ_p& out, size_t n_read, const uint8_t* buf) {
    static ZZ tmp(INIT_SIZE, _256_BITS); // Allocate once
    ZZFromBytes(tmp, buf, n_read);
    conv(out, tmp);
}

void buff_to_mcl(mcl::Fr& out, size_t n_read, const uint8_t* buf) {
    out.setLittleEndianMod(buf, n_read);
}

/******************************************************************************/
/**************************** MPZ <--> NTL ************************************/
/******************************************************************************/
void mpz_to_ZZ_p(ZZ_p& out, const mpz_class& in) {

    size_t n_wrote = 42; // n_wrote should never be 42;
    mpz_to_buff(BUF, &n_wrote, in);

    char* mpz_str = print_to_string_mpz1(in);
    std::cout << "mpz_to_ZZ_p is serialising the mpz_t:" << mpz_str << "\n";
    assert(n_wrote != 42);

    buff_to_ntl(out, n_wrote, BUF);

    // clear buff
    memset((void*)BUF, 0, N_BYTES_256_BITS);
}

ZZ_p mpz_to_new_ZZ_p(const mpz_class& in) {
    ZZ_p out(INIT_ALLOC);
    mpz_to_ZZ_p(out, in);
    return out;
}
void ZZ_p_to_mpz(mpz_class& out, const ZZ_p& in) {
    ZZ inzz; inzz = rep(in); // TODO is conv<ZZ>(in) better?
    long n_wrote = 42;

    ntl_to_buff(BUF, &n_wrote, inzz);
    assert(n_wrote >= 0 && n_wrote <= 32); // 256 bits (2^8) == 2^5 Bytes

    buff_to_mpz(out, n_wrote, BUF);
    memset((void*)BUF, 0, N_BYTES_256_BITS);
}

/******************************************************************************/
/**************************** NTL |--> MCL ************************************/
/******************************************************************************/

// a version that's better if all Fr's are allocated ahead of time
void ZZ_p_to_Fr(mcl::Fr& out, const ZZ_p& in) {
    ZZ inzz; inzz = rep(in); // TODO is conv<ZZ>(in) better?
    long n_wrote = 42;
    ntl_to_buff(BUF, &n_wrote, inzz);
    buff_to_mcl(out, n_wrote, BUF);
    memset((void*)BUF, 0, N_BYTES_256_BITS);
}

mcl::Fr ZZ_p_to_new_Fr(const ZZ_p& in) {
    mcl::Fr out;
    ZZ_p_to_Fr(out, in);
    return out;
}

/******************************************************************************/
/**************************** GMP |--> MCL ************************************/
/******************************************************************************/

// a version that's better if all Fr's are allocated ahead of time
void mpz_to_Fr(mcl::Fr& out, const mpz_class& in) {
    size_t n_wrote = 42; // n_wrote should never be 42;
    mpz_to_buff(BUF, &n_wrote, in);
    assert(n_wrote <= 32); // 256 bits (2^8) == 2^5 Bytes

    buff_to_mcl(out, n_wrote, BUF);
    memset((void*)BUF, 0, N_BYTES_256_BITS);
}

mcl::Fr mpz_to_new_Fr(const mpz_class& in) {
    mcl::Fr out;
    mpz_to_Fr(out, in);
    return out;
}

