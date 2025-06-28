#include <iostream>
#include <cstdlib>
#include <vector>
#include <cassert>
#include <cmath>

typedef __int128_t i128;
using namespace std;

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

vector<i128>
round_vec(vector<double> v) {
    vector<i128> v_rounded(v.size());
    for (size_t i = 0; i < v.size(); i++)
        v_rounded.at(i) = static_cast<i128>(round(v.at(i)));
    return v_rounded;
}

vector<i128>
add_vec(vector<i128> v1, vector<i128> v2, i128 q) {
    vector<i128> v(v1.size());
    for (size_t i = 0; i < v1.size(); i++) {
        assert(v.at(i) == 0);
        v.at(i) = v1.at(i) + v2.at(i);
    }
    return v;
}

void
test_round_vec() {
    vector<double> v{-3.4, -5.6, 2.3, 7.7};
    vector<i128> v_rounded = round_vec(v);
    cout << "v: ";
    for (auto el : v_rounded) {
        print_i128(el);
        cout << ", ";
    }
    cout << "\n";
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
    test_round_vec();
    return EXIT_SUCCESS;
}