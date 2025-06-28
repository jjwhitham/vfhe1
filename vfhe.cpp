#include <iostream>
#include <cstdlib>
#include <vector>
#include <cassert>

typedef __int128_t i128;
using namespace std;

vector<vector<i128>>
scalar_mat_mult(i128 s, vector<vector<i128>> M, i128 q) {
    // """Multiplies each element of matrix M by scalar s (mod q)."""
    // return [[int(s * elem) % q for elem in row] for row in M]
    for (size_t i = 0; i < M.size(); i++) {
        for (size_t j = 0; j < M.at(0).size(); j++) {
            M.at(i).at(j) = M.at(i).at(j) * s % q;
        }
    }
    return M;
}

vector<i128>
scalar_vec_mult(i128 s, vector<i128> v, i128 q) {
    // """Multiplies each element of matrix M by scalar s (mod q)."""
    // return [[int(s * elem) % q for elem in row] for row in M]
    for (size_t i = 0; i < v.size(); i++) {
        v.at(i) = v.at(i) * s % q;
    }
    return v;
}

void print_i128(i128 n) {
    if (n == 0) {
        cout << "0";
        return;
    }
    bool neg = false;
    if (n < 0) {
        neg = true;
        n = -n;
    }
    char buf[50];
    int i = 49;
    buf[i--] = '\0';
    while (n > 0) {
        buf[i--] = '0' + (n % 10);
        n /= 10;
    }
    if (neg) buf[i--] = '-';
    cout << (buf + i + 1);
}

void test_scalar_mat_mult() {
    int rows = 3;
    vector<vector<i128>> M(rows);
    int cols = 2;
    for (int i = 0; i < rows; i++) {
        vector<i128> row(cols);
        M.at(i) = row;
    }
    int k = 1;
    for (int i = 0; i < rows; i++)
        for (int j = 0; j < cols; j++)
            M.at(i).at(j) = k++;
    int s = 2;
    M = scalar_mat_mult(s, M, 13);
    cout << "M[2][1]: ";
    print_i128(M.at(2).at(1));
    cout << "\n";
    assert(M.at(2).at(1) == 12);
}

int main() {
    test_scalar_mat_mult();
    return EXIT_SUCCESS;
}