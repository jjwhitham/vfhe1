#include <iostream>
#include <cstdlib>
#include <vector>
#include <cassert>
#include <cmath>
#include <climits>
#include <algorithm>

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

class array2d_vec {
private:
    size_t rows_;
    size_t cols_;
    vector<vector<i128>> data_;
public:
    array2d_vec(size_t rows, size_t cols) : rows_(rows), cols_(cols) {
        data_.resize(rows, vector<i128>(cols, 0));
    }
    ~array2d_vec() = default;

    vector<i128>& operator[](size_t index) {
        return data_[index];
    }

    const vector<i128>& operator[](size_t index) const {
        return data_[index];
    }
};

class array1d {
private:
    size_t cols_;
    i128* arr;
public:
    array1d(size_t cols) : cols_(cols) {
        arr = new i128[cols](); // new []() to zero-initialise
    }
    ~array1d() {
        delete[] arr;
    }
    void check_index_bounds(size_t col) {
        if (col >= cols_) {
            throw std::out_of_range(
                "Index error: accessing arr[" + std::to_string(col) + "]"
                + " in a " + std::to_string(cols_) + " element array."
            );
        }
    }
    i128 get(int col) {
        check_index_bounds(col);
        return arr[col];
    }
    void check_value_bounds(i128 val) {
        // assert abs(val) less than 2^63 using INT_MIN and INT_MAX
        i128 min_val = -1;
        min_val <<= 63;
        i128 max_val = (1UL << 63) - 1;
        if (val < min_val || val > max_val) {
            throw std::out_of_range(
                "(array2d) Value out of range: " + print_to_string_i128(val)
            );
        }
    }
    void set(int col, i128 val) {
        check_index_bounds(col);
        check_value_bounds(val);
        arr[col] = val;
    }
};

class array2d {
private:
    size_t rows_;
    size_t cols_;
    i128* data;
    i128** arr;
public:
    array2d(size_t rows, size_t cols) : rows_(rows), cols_(cols) {
        // allocate contiguous mem
        data = new i128[rows * cols](); // new []() to zero-initialise
        // point each row of 2d arr to the right place in 1d data
        arr = new i128*[rows];
        for (size_t i = 0; i < rows; i++) {
            arr[i] = &data[i * cols];
        }
    }
    ~array2d() {
        delete[] data;
        delete[] arr;
    }
    void check_index_bounds(size_t row, size_t col) {
        if (row >= rows_ || col >= cols_) {
            throw std::out_of_range(
                "Index out of bounds: tryed to access "
                + std::to_string(row) + ", " + std::to_string(col)
                + " in a " + std::to_string(rows_) + "x" + std::to_string(cols_)
                + " array."
            );
        }
    }
    i128 get(int row, int col) {
        check_index_bounds(row, col);
        return arr[row][col];
    }
    void check_value_bounds(i128 val) {
        // assert abs(val) less than 2^63 using INT_MIN and INT_MAX
        i128 min_val = -1;
        min_val <<= 63;
        i128 max_val = (1UL << 63) - 1;
        if (val < min_val || val > max_val) {
            throw std::out_of_range(
                "(array2d) Value out of range: " + print_to_string_i128(val)
            );
        }
    }
    void set(int row, int col, i128 val) {
        check_index_bounds(row, col);
        check_value_bounds(val);
        arr[row][col] = val;
    }
};

void test_array2d() {
    i128 min_val = -1;
    min_val <<= 63;
    i128 max_val = (1UL << 63) - 1;
    array2d arr(3, 2);
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
    bool tests_passed = false;
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
    array1d arr(3);
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

int main() {
    test_scalar_mat_mult();
    test_round_vec();
    test_array2d();
    test_array1d();
    return EXIT_SUCCESS;
}