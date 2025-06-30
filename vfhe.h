#pragma once

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

template<typename T>
class array1d {
private:
    size_t size_;
    T* arr;
    // Proxy class for operator[]
    class Proxy {
        array1d<T>& parent;
        size_t idx;
    public:
        Proxy(array1d<T>& parent_, size_t idx_);
        // Assignment operator with value check
        Proxy& operator=(T val);
        // Conversion operator for reading value
        operator T() const;
    };
public:
    array1d(size_t size);
    ~array1d();
    void check_index_bounds(size_t n) const;
    T get(size_t n);
    // Use Proxy for operator[]
    Proxy operator[](size_t idx);
    // For const access
    T operator[](size_t idx) const;
    void check_value_bounds(T val);
    void set(int n, T val);
    size_t size() const;
};


class array2d {
private:
    size_t rows_;
    size_t cols_;
    i128* data;
    i128** arr;
    // Proxy for a row, to enable arr[row][col] with checks
    class RowProxy {
        array2d& parent;
        size_t row;
    public:
        RowProxy(array2d& parent_, size_t row_);
        // Assignment and access with value/index checks
        class ColProxy {
            array2d& parent;
            size_t row, col;
        public:
            ColProxy(array2d& parent_, size_t row_, size_t col_);
            // Assignment operator with checks
            ColProxy& operator=(i128 val);
            // Conversion operator for reading value
            operator i128() const;
        };
        ColProxy operator[](size_t col);
        // For const access
        i128 operator[](size_t col) const;
    };

public:
    array2d(size_t rows, size_t cols) ;
    ~array2d();
    void check_index_bounds(size_t row, size_t col) const;
    i128 get(size_t row, size_t col) const;
    void check_value_bounds(i128 val) const;
    void set(size_t row, size_t col, i128 val);
    size_t size() const;
    size_t size(size_t row) const;
    // Proxy for arr[row][col] with checks
    RowProxy operator[](size_t row);
    // For const access
    const RowProxy operator[](size_t row) const;
};


void print_i128(i128 n);

vector<vector<i128>>
scalar_mat_mult(i128 s, vector<vector<i128>> M, i128 q);

array2d
scalar_mat_mult(i128 s, array2d M, i128 q);

vector<i128>
scalar_vec_mult(i128 s, vector<i128> v, i128 q);

array1d<i128>
scalar_vec_mult(i128 s, array1d<i128> v, i128 q);

vector<i128>
round_vec(vector<double> v);

vector<i128>
add_vec(vector<i128> v1, vector<i128> v2, i128 q);

// class array2d_vec {
// private:
//     size_t rows_;
//     size_t cols_;
//     vector<vector<i128>> data_;
// public:
//     array2d_vec(size_t rows, size_t cols);
//     ~array2d_vec();

//     vector<i128>& operator[](size_t index);

//     const vector<i128>& operator[](size_t index) const;
// };
