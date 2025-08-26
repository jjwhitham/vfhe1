#include <iostream>
#include <gmpxx.h>
#include "gmp.h"

int main() {
    mpz_class m = 5;
    std::cout << "\nModulus m: ";
    char buf[100];
    buf[99] = '\0';
    gmp_sprintf(buf, "%Zd", m.get_mpz_t());
    std::cout << buf << "\n\n";
}