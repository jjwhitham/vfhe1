#include <iostream>
#include <cstdlib>
#include <vector>
#include <cassert>
#include <cmath>
#include <climits>
#include <algorithm>

typedef __int128_t i128;
using namespace std;

string print_to_string_i128(i128 n);

class array1d {
private:
    size_t size_;
    i128* arr;
    // Proxy class for operator[]
    class Proxy {
        array1d& parent;
        size_t idx;
    public:
        Proxy(array1d& parent_, size_t idx_);
        // Assignment operator with value check
        Proxy& operator=(i128 val);
        // Conversion operator for reading value
        operator i128() const;
    };
public:
    array1d(size_t size);
    ~array1d();
    void check_index_bounds(size_t n) const;
    i128 get(size_t n);
    // Use Proxy for operator[]
    Proxy operator[](size_t idx);
    // For const access
    i128 operator[](size_t idx) const;
    void check_value_bounds(i128 val);
    void set(int n, i128 val);
    size_t size() const;
};

class array2d {
private:
    size_t rows_;
    size_t cols_;
    i128* data;
    i128** arr;
public:
    array2d(size_t rows, size_t cols);
    ~array2d();
    void check_index_bounds(size_t row, size_t col);
    i128 get(int row, int col);
    void check_value_bounds(i128 val);
    void set(int row, int col, i128 val);
    size_t size();
    size_t size(size_t row);
};

void print_i128(i128 n);

vector<vector<i128>>
scalar_mat_mult(i128 s, vector<vector<i128>> M, i128 q);

array2d
scalar_mat_mult(i128 s, array2d M, i128 q);

vector<i128>
scalar_vec_mult(i128 s, vector<i128> v, i128 q);

array1d
scalar_vec_mult(i128 s, array1d v, i128 q);

vector<i128>
round_vec(vector<double> v);

vector<i128>
add_vec(vector<i128> v1, vector<i128> v2, i128 q);

class array2d_vec {
private:
    size_t rows_;
    size_t cols_;
    vector<vector<i128>> data_;
public:
    array2d_vec(size_t rows, size_t cols);
    ~array2d_vec();

    vector<i128>& operator[](size_t index);

    const vector<i128>& operator[](size_t index) const;
};
