// TODO array1d/2d<double>
#include <iostream>
#include <cstdlib>
#include <vector>
#include <cassert>
#include <cmath>
#include <climits>
#include <algorithm>
#include "vfhe.h"

// typedef __int128_t i128;
// using namespace std;

string print_to_string_i128(i128 n) {
    if (n == 0) {
        return "0";
    }
    bool neg = false;
    if (n < 0) {
        neg = true;
        n = -n;
    }
    string buf;
    while (n > 0) {
        buf += '0' + (n % 10);
        n /= 10;
    }
    if (neg) buf += '-';
    reverse(buf.begin(), buf.end());
    return buf;
}
// template class array1d<i128>;
template<typename T>
array1d<T>::array1d(size_t size) : size_(size) {
    arr = new T[size]();
}
template<typename T>
array1d<T>::~array1d() {
    delete[] arr;
}
template<typename T>
void array1d<T>::check_index_bounds(size_t n) const {
    if (n >= size_) {
        throw std::out_of_range(
            "Index error: accessing arr[" + std::to_string(n) + "]"
            + " in a " + std::to_string(size_) + " element array."
        );
    }
}
template<typename T>
T array1d<T>::get(size_t n) {
    check_index_bounds(n);
    return arr[n];
}
// Use Proxy for operator[]
template<typename T>
typename array1d<T>::Proxy array1d<T>::operator[](size_t idx) {
    return Proxy(*this, idx);
}
// For const access
template<typename T>
T array1d<T>::operator[](size_t idx) const {
    check_index_bounds(idx);
    return arr[idx];
}
template<>
void array1d<i128>::check_value_bounds(i128 val) {
    i128 min_val = -1;
    min_val <<= 63;
    i128 max_val = (1UL << 63) - 1;
    if (val < min_val || val > max_val) {
        throw std::out_of_range(
            "(array2d) Value out of range: " + print_to_string_i128(val)
        );
    }
}
template<typename T>
void array1d<T>::check_value_bounds(T val) {
    cout << "check_value_bounds not implemented for type T\n";
}
template<typename T>
void array1d<T>::set(int n, T val) {
    check_index_bounds(n);
    check_value_bounds(val);
    arr[n] = val;
}
template<typename T>
size_t array1d<T>::size() const {
    return size_;
}

template<typename T>
array1d<T>::Proxy::Proxy(array1d<T>& parent_, size_t idx_) : parent(parent_), idx(idx_) {}
// Assignment operator with value check
template<typename T>
typename array1d<T>::Proxy& array1d<T>::Proxy::operator=(T val) {
    parent.check_index_bounds(idx);
    parent.check_value_bounds(val);
    parent.arr[idx] = val;
    return *this;
}
// Conversion operator for reading value
template<typename T>
array1d<T>::Proxy::operator T() const {
    parent.check_index_bounds(idx);
    return parent.arr[idx];
}


array2d::RowProxy::RowProxy(array2d& parent_, size_t row_) : parent(parent_), row(row_) {}

// Assignment and access with value/index checks
array2d::RowProxy::ColProxy::ColProxy(array2d& parent_, size_t row_, size_t col_)
        : parent(parent_), row(row_), col(col_) {}

    // Assignment operator with checks
   array2d::RowProxy::ColProxy& array2d::RowProxy::ColProxy::operator=(i128 val) {
        parent.check_index_bounds(row, col);
        parent.check_value_bounds(val);
        parent.arr[row][col] = val;
        return *this;
    }

// Conversion operator for reading value
array2d::RowProxy::ColProxy::operator i128() const {
    parent.check_index_bounds(row, col);
    return parent.arr[row][col];
}

array2d::RowProxy::ColProxy array2d::RowProxy::operator[](size_t col) {
    return ColProxy(parent, row, col);
}
// For const access
i128 array2d::RowProxy::operator[](size_t col) const {
    parent.check_index_bounds(row, col);
    return parent.arr[row][col];
}

array2d::array2d(size_t rows, size_t cols) : rows_(rows), cols_(cols) {
    data = new i128[rows * cols]();
    arr = new i128*[rows];
    for (size_t i = 0; i < rows; i++) {
        arr[i] = &data[i * cols];
    }
}
array2d::~array2d() {
    delete[] data;
    delete[] arr;
}
void array2d::check_index_bounds(size_t row, size_t col) const {
    if (row >= rows_ || col >= cols_) {
        throw std::out_of_range(
            "Index out of bounds: tried to access "
            + std::to_string(row) + ", " + std::to_string(col)
            + " in a " + std::to_string(rows_) + "x" + std::to_string(cols_)
            + " array."
        );
    }
}
i128 array2d::get(size_t row, size_t col) const {
    check_index_bounds(row, col);
    return arr[row][col];
}
void array2d::check_value_bounds(i128 val) const {
    i128 min_val = -1;
    min_val <<= 63;
    i128 max_val = (1UL << 63) - 1;
    if (val < min_val || val > max_val) {
        throw std::out_of_range(
            "(array2d) Value out of range: " + print_to_string_i128(val)
        );
    }
}
void array2d::set(size_t row, size_t col, i128 val) {
    check_index_bounds(row, col);
    check_value_bounds(val);
    arr[row][col] = val;
}
size_t array2d::size() const {
    return rows_;
}
size_t array2d::size(size_t row) const {
    check_index_bounds(row, 0);
    return cols_;
}

// Proxy for arr[row][col] with checks
array2d::RowProxy array2d::operator[](size_t row) {
    return RowProxy(*this, row);
}
// For const access
const array2d::RowProxy array2d::operator[](size_t row) const {
    return RowProxy(const_cast<array2d&>(*this), row);
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

array2d
scalar_mat_mult(i128 s, array2d M, i128 q) {
    // """Multiplies each element of matrix M by scalar s (mod q)."""
    // return [[int(s * elem) % q for elem in row] for row in M]
    for (size_t i = 0; i < M.size(); i++) {
        for (size_t j = 0; j < M.size(0); j++) {
            i128 new_val = M.get(i, j) * s % q;
            M.set(i, j, new_val);
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

array1d<i128>
scalar_vec_mult(i128 s, array1d<i128> v, i128 q) {
    // """Multiplies each element of matrix M by scalar s (mod q)."""
    // return [[int(s * elem) % q for elem in row] for row in M]
    for (size_t i = 0; i < v.size(); i++) {
        i128 new_val = v.get(i) * s % q;
        v.set(i, new_val);
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
        v.at(i) = (v1.at(i) + v2.at(i)) % q;
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

// class array2d_vec {
// private:
//     size_t rows_;
//     size_t cols_;
//     vector<vector<i128>> data_;
// public:
//     array2d_vec(size_t rows, size_t cols) : rows_(rows), cols_(cols) {
//         data_.resize(rows, vector<i128>(cols, 0));
//     }
//     ~array2d_vec() = default;

//     vector<i128>& operator[](size_t index) {
//         return data_[index];
//     }

//     const vector<i128>& operator[](size_t index) const {
//         return data_[index];
//     }
// };

void test_array2d() {
    i128 min_val = -1;
    min_val <<= 63;
    i128 max_val = (1UL << 63) - 1;
    array2d arr(3, 2);
    // check values are zero initialised
    for (size_t i = 0; i < arr.size(); i++) {
        for (size_t j = 0; j < arr.size(i); j++) {
            assert(arr.get(i, j) == 0);
        }
    }
    // set some values
    arr.set(0, 0, 1);
    arr.set(1, 1, 2);
    arr.set(2, 0, 3);
    cout << "arr[0][0]: ";
    print_i128(arr.get(0, 0));
    cout << "\n";
    cout << "arr[1][1]: ";
    print_i128(arr.get(1, 1));
    cout << "\n";
    cout << "arr[2][0]: ";
    print_i128(arr.get(2, 0));
    cout << "\n";
    assert(arr.get(0, 0) == 1);
    assert(arr.get(1, 1) == 2);
    assert(arr.get(2, 0) == 3);
    // Test out of bounds
    // bool tests_passed = false;
    try {
        arr.get(3, 0);
    } catch (const std::out_of_range& e) {
        cout << "Caught expected exception: " << e.what() << "\n";
    }
    try {
        arr.get(0, 2);
    } catch (const std::out_of_range& e) {
        cout << "Caught expected exception: " << e.what() << "\n";
    }
    try {
        arr.set(0, 2, 5);
    } catch (const std::out_of_range& e) {
        cout << "Caught expected exception: " << e.what() << "\n";
    }
    try {
        arr.set(3, 0, 5);
    } catch (const std::out_of_range& e) {
        cout << "Caught expected exception: " << e.what() << "\n";
    }
    try {
        arr.set(0, 0, min_val - 1);
    } catch (const std::out_of_range& e) {
        cout << "Caught expected exception: " << e.what() << "\n";
    }
    try {
        arr.set(0, 0, max_val + 1);
    } catch (const std::out_of_range& e) {
        cout << "Caught expected exception: " << e.what() << "\n";
    }
    try {
        arr.set(0, 0, max_val); // should not throw
    } catch (const std::out_of_range& e) {
        cout << "Caught unexpected exception: " << e.what() << "\n";
    }
    try {
        arr.set(0, 0, min_val); // should not throw
    } catch (const std::out_of_range& e) {
        cout << "Caught unexpected exception: " << e.what() << "\n";
    }
    // if (!tests_passed)
    //     throw runtime_error("test_array2d failed");
    // else
    cout << "All tests passed for array2d.\n";
};

void test_array1d() {
    i128 min_val = -1;
    min_val <<= 63;
    i128 max_val = (1UL << 63) - 1;

    array1d<i128> arr(3);
    // check values are zero initialised
    for (size_t i = 0; i < arr.size(); i++) {
        assert(arr.get(i) == 0);
    }
    arr.set(0, 1);
    arr.set(1, 2);
    arr.set(2, 3);
    cout << "arr[0]: ";
    print_i128(arr.get(0));
    cout << "\n";
    cout << "arr[1]: ";
    print_i128(arr.get(1));
    cout << "\n";
    cout << "arr[2]: ";
    print_i128(arr.get(2));
    cout << "\n";
    assert(arr.get(0) == 1);
    assert(arr.get(1) == 2);
    assert(arr.get(2) == 3);
    try {
        arr.get(3);
    } catch (const std::out_of_range& e) {
        cout << "Caught expected exception: " << e.what() << "\n";
    }
    try {
        arr.set(3, 5);
    } catch (const std::out_of_range& e) {
        cout << "Caught expected exception: " << e.what() << "\n";
    }
    try {
        arr.set(0, min_val - 1);
    } catch (const std::out_of_range& e) {
        cout << "Caught expected exception: " << e.what() << "\n";
    }
    try {
        arr.set(0, max_val + 1);
    } catch (const std::out_of_range& e) {
        cout << "Caught expected exception: " << e.what() << "\n";
    }
    try {
        arr.set(0, max_val); // should not throw
    } catch (const std::out_of_range& e) {
        cout << "Caught unexpected exception: " << e.what() << "\n";
    }
    try {
        arr.set(0, min_val); // should not throw
    } catch (const std::out_of_range& e) {
        cout << "Caught unexpected exception: " << e.what() << "\n";
    }
    cout << "All tests passed for array1d.\n";
}

void test_array1d_indexing() {
    array1d<i128> arr(5);
    for (size_t i = 0; i < arr.size(); i++) {
        arr[i] = i + 1; // Using operator[] for assignment
    }
    for (size_t i = 0; i < arr.size(); i++) {
        assert(arr[i] == i + 1); // Using operator[] for access
    }
    cout << "All tests passed for array1d indexing.\n";
}

void test_array2d_indexing() {
    array2d arr(3, 3);
    for (size_t i = 0; i < arr.size(); i++) {
        for (size_t j = 0; j < arr.size(i); j++) {
            arr[i][j] = i * 3 + j + 1; // Using operator[] for assignment
        }
    }
    for (size_t i = 0; i < arr.size(); i++) {
        for (size_t j = 0; j < arr.size(i); j++) {
            assert(arr[i][j] == i * 3 + j + 1); // Using operator[] for access
        }
    }
    cout << "All tests passed for array2d indexing.\n";
}

// int main() {
//     test_scalar_mat_mult();
//     test_round_vec();
//     test_array1d();
//     test_array2d();
//     test_array1d_indexing();
//     test_array2d_indexing();
//     return EXIT_SUCCESS;
// }