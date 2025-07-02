/*
Setup:
rF, rG, sH
- veri_vec * rgsw_vec
    - scalar * rgsw

Verfication:
rx' = rFx + rGx
su = sHu
*/

// TODO array1d/2d<double>
#include <iostream>
#include <cstdlib>
#include <vector>
#include <cassert>
#include <cmath>
#include <climits>
#include <algorithm>

typedef __int128_t i128;
using namespace std;
// global variable for modulus
constexpr i128 GROUP_MODULUS = 11;
constexpr i128 GENERATOR = 6;

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

template<typename T, typename Derived>
class array1d {
private:
    size_t size_;
    T* arr;
public:
    array1d(size_t size) : size_(size) {
        arr = new T[size]();
    }
    ~array1d() = default;
    //  {

        // delete[] arr;
    // }
    void check_index_bounds(size_t n) const {
        if (n >= size_) {
            throw std::out_of_range(
                "Index error: accessing arr[" + std::to_string(n) + "]"
                + " in a " + std::to_string(size_) + " element array."
            );
        }
    }
    T& get(size_t n) const {
        check_index_bounds(n);
        return arr[n];
    }
    void check_value_bounds(T val) {
        if constexpr (std::is_same_v<T, i128>) {
            i128 min_val = -1;
            min_val <<= 63;
            i128 max_val = (1UL << 63) - 1;
            if (val < min_val || val > max_val) {
                throw std::out_of_range(
                    "(array2d) Value out of range: " + print_to_string_i128(val)
                );
            }
            // else
                // cout << "check_value_bounds not implemented for type T\n";
        }
    }
    void set(int n, T val) {
        check_index_bounds(n);
        check_value_bounds(val);
        arr[n] = val;
    }
    size_t size() const {
        return size_;
    }
    // when derived class is used, this is called recursively
    Derived operator*(const Derived& other) const {
        size_t N = size();
        Derived res(N);
        for (size_t i = 0; i < N; i++) {
            auto val = (get(i) * other.get(i));
            if constexpr (std::is_same_v<T, i128>)
                val %= GROUP_MODULUS;
            res.set(i, val);
        }
        return res;
    }
    // Scalar multiplication, for:
    // rlwe_decom_vec * rgsw_mat
    Derived operator*(const i128& scalar) const {
        size_t N = size();
        Derived res(N);
        for (size_t i = 0; i < N; i++) {
            // multiply each element by scalar
            auto val = (get(i) * scalar) % GROUP_MODULUS;
            res.set(i, val);
        }
        return res;
    }
    // when derived class is used, this is called recursively
    Derived operator+(const Derived& other) const {
        size_t N = size();
        Derived res(N);
        for (size_t i = 0; i < N; i++) {
            auto val = get(i) + other.get(i);
            res.set(i, val);
        }
        return res;
    }
    Derived pow1(const Derived& other) const {
        size_t N = size();
        Derived res(N);
        if constexpr (std::is_same_v<T, i128>) {
            for (size_t i = 0; i < N; i++) {
                auto val = pow1_(get(i), other.get(i));
                res.set(i, val);
            }
            return res;
        } else {
            for (size_t i = 0; i < N; i++) {
                auto val = get(i).pow1(other.get(i));
                res.set(i, val);
            }
            return res;
        }
    }
    i128 pow1_(i128 self, i128 other) const {
        i128 square = GENERATOR;
        while (other != 0) {
            if ((other & 1) == 1)
                self = (self * square) % GROUP_MODULUS;
            square = (square * square) % GROUP_MODULUS;
            other >>= 1;
        }
        return self;
    }

    T* begin() { return arr; }
    T* end() { return arr + size_; }
    const T* begin() const { return arr; }
    const T* end() const { return arr + size_; }
};

template<typename T>
class array2d {
private:
    size_t rows_;
    size_t cols_;
    T* data;
    T** arr;
public:
    array2d(size_t rows, size_t cols) : rows_(rows), cols_(cols) {
        data = new T[rows * cols]();
        arr = new T*[rows];
        for (size_t i = 0; i < rows; i++) {
            arr[i] = &data[i * cols];
        }
    }
    ~array2d() {
        delete[] data;
        delete[] arr;
    }
    void check_index_bounds(size_t row, size_t col) const {
        if (row >= rows_ || col >= cols_) {
            throw std::out_of_range(
                "Index out of bounds: tried to access "
                + std::to_string(row) + ", " + std::to_string(col)
                + " in a " + std::to_string(rows_) + "x" + std::to_string(cols_)
                + " array."
            );
        }
    }
    T& get(size_t row, size_t col) const {
        check_index_bounds(row, col);
        return arr[row][col];
    }
    void check_value_bounds(i128 val) const {
        i128 min_val = -1;
        min_val <<= 63;
        i128 max_val = (1UL << 63) - 1;
        if (val < min_val || val > max_val) {
            throw std::out_of_range(
                "(array2d) Value out of range: " + print_to_string_i128(val)
            );
        }
    }
    void check_value_bounds(T val) {
        cout << "check_value_bounds not implemented for type T\n";
    }
    void set(size_t row, size_t col, T val) {
        check_index_bounds(row, col);
        check_value_bounds(val);
        arr[row][col] = val;
    }
    size_t size() const {
        return rows_;
    }
    size_t size(size_t row) const {
        check_index_bounds(row, 0);
        return cols_;
    }
};

class poly : public array1d<i128, poly> {
private:
public:
    static constexpr size_t DEFAULT_N = 4;
    poly(size_t N) : array1d<i128, poly>(N) {
        for (size_t i = 0; i < N; i++) {
            set(i, i + 1); // initialize coefficients to zero
        }
    }
    poly() : array1d<i128, poly>(DEFAULT_N) {
        for (size_t i = 0; i < DEFAULT_N; i++) {
            set(i, i + 1); // initialize coefficients to zero
        }
    }
    auto& get_coeff(size_t n) const {
        return get(n);
    }
    // poly operator*(const poly& other) const {
    //     size_t N = size();
    //     poly res(N);
    //     for (size_t i = 0; i < N; i++) {
    //         i128 val = (*this).get(i) * other.get(i);
    //         res.set(i, val);
    //     }
    //     return res;
    // }
};

/* a tuple of polys */
class rlwe : public array1d<poly, rlwe> {
private:
public:
    static constexpr size_t DEFAULT_N = 2;
    rlwe(size_t N) : array1d<poly, rlwe>(N) {}
    rlwe() : array1d<poly, rlwe>(DEFAULT_N) {}
    auto& get_poly(size_t n) const {
        return get(n);
    }
};

class rlwe_vec : public array1d<rlwe, rlwe_vec> {
private:
public:
    static constexpr size_t DEFAULT_N = 4;
    rlwe_vec(size_t N) : array1d<rlwe, rlwe_vec>(N) {}
    rlwe_vec() : array1d<rlwe, rlwe_vec>(DEFAULT_N) {}
    auto& get_rlwe(size_t n) const {
        return get(n);
    }
};

class rlwe_decomp : public array1d<poly, rlwe_decomp> {
private:
public:
    static constexpr size_t DEFAULT_N = 4;
    rlwe_decomp(size_t N) : array1d<poly, rlwe_decomp>(N) {}
    rlwe_decomp() : array1d<poly, rlwe_decomp>(DEFAULT_N) {}
    ~rlwe_decomp() {}
};

class rlwe_decomp_vec : public array1d<rlwe_decomp, rlwe_decomp_vec> {
private:
public:
    static constexpr size_t DEFAULT_N = 3;
    rlwe_decomp_vec(size_t N) : array1d<rlwe_decomp, rlwe_decomp_vec>(N) {}
    rlwe_decomp_vec() : array1d<rlwe_decomp, rlwe_decomp_vec>(DEFAULT_N) {}
    ~rlwe_decomp_vec() {}
};

// class veri_vec_scalar : public array1d<i128, veri_vec_scalar> {
// private:
// public:
//     // TODO mult for veri_vec_scalar * rgsw_mat
//     rlwe operator*(const rlwe_decomp& other) const {
//         size_t N = size();
//         assert(N == other.size());
//         rlwe res;
//         poly p0 = res.get(0);
//         poly p1 = res.get(1);
//         for (size_t i = 0; i < N; i++) {
//             auto val0 = get(i).get(0) * other.get(i);
//             auto val1 = get(i).get(1) * other.get(i);
//             p0 = p0 + val0;
//             p1 = p1 + val1;
//         }
//         res.set(0, p0);
//         res.set(1, p1);
//         return res;
//     }

    // TODO pow func for g^rFx, etc
    // rlwe_vec pow(const rgsw_mat& other) const {
    //     size_t rows = other.size();
    //     size_t cols = other.size(0);
    //     assert(rows == size());
    //     rlwe_vec res;
    //     rlwe sum;
    //     for (size_t j = 0; j < cols; j++) {
    //         // for each element in the row
    //         for (size_t i = 0; i < rows; i++) {
    //             auto val = get(i).pow(other.get(i, j)).pow(other.get(j));
    //             sum = sum * val;
    //         }
    //         res.set(i, sum);
    //     }
    //     return res;
    // }
// };

class rgsw : public array1d<rlwe, rgsw> {
private:
public:
    static constexpr size_t DEFAULT_N = 4;
    rgsw(size_t N) : array1d<rlwe, rgsw>(N) {}
    rgsw() : array1d<rlwe, rgsw>(DEFAULT_N) {}
    ~rgsw() {}

    rlwe operator*(const rlwe_decomp& other) const {
        size_t N = size();
        assert(N == other.size());
        rlwe res;
        poly p0 = res.get(0);
        poly p1 = res.get(1);
        for (size_t i = 0; i < N; i++) {
            auto val0 = get(i).get(0) * other.get(i);
            auto val1 = get(i).get(1) * other.get(i);
            p0 = p0 + val0;
            p1 = p1 + val1;
        }
        res.set(0, p0);
        res.set(1, p1);
        return res;
    }
    rlwe pow(const rlwe_decomp& other) const {
        size_t N = size();
        assert(N == other.size());
        rlwe res;
        poly p0 = res.get(0);
        poly p1 = res.get(1);
        for (size_t i = 0; i < N; i++) {
            auto val0 = get(i).get(0).pow1(other.get(i));
            auto val1 = get(i).get(1).pow1(other.get(i));
            p0 = p0 * val0;
            p1 = p1 * val1;
        }
        res.set(0, p0);
        res.set(1, p1);
        return res;
    }
};

class rgsw_mat : public array2d<rgsw> {
private:
public:
    rgsw_mat() : array2d<rgsw>(3, 3) {}
    ~rgsw_mat() {}
    rlwe_vec operator*(const rlwe_decomp_vec& other) const {
        size_t rows = size();
        size_t cols = size(0);
        assert(cols == other.size());
        rlwe_vec res;
        rlwe sum;
        for (size_t i = 0; i < rows; i++) {
            // for each element in the row
            for (size_t j = 0; j < cols; j++) {
                auto val = get(i, j) * other.get(j);
                sum = sum + val;
            }
            res.set(i, sum);
        }
        return res;
    }
    rlwe_vec pow(const rlwe_decomp_vec& other) const {
        size_t rows = size();
        size_t cols = size(0);
        assert(cols == other.size());
        rlwe_vec res;
        rlwe sum;
        for (size_t i = 0; i < rows; i++) {
            // for each element in the row
            for (size_t j = 0; j < cols; j++) {
                auto val = get(i, j).pow(other.get(j));
                sum = sum * val;
            }
            res.set(i, sum);
        }
        return res;
    }
};

int main() {
    rlwe r;
    r.get(0).set(0, 21);
    r.get(1).set(0, 2);

    poly p = r.get(0) * r.get(1);
    // for (int i = 0; i < 4; i++) {
    //     cout << print_to_string_i128(p.get(i));
    //     cout << ", ";
    // }
    // cout << "\n";


    rlwe_vec rv;
    poly p1 = rv.get(0).get(0);
    p1.set(0, 21);
    poly p2 = rv.get(0).get(1);
    p2.set(0, 2);
    poly p3 = p1 * p2;
    // for (int i = 0; i < 4; i++) {
    //     cout << print_to_string_i128(p3.get(i));
    //     cout << ", ";
    // }
    // cout << "\n";
    // for (int i = 0; i < 4; i++) {
    //     cout << print_to_string_i128(p3.get(i));
    //     cout << ", ";
    // }
    // cout << "\n";
    rv.get(1).set(0, p3);
    rlwe_vec r1;
    rlwe_vec r2;
    rlwe_vec r3 = r1 * r2;

    rgsw_mat F;
    rlwe_decomp_vec x;
    auto Fx1 = F * x;
    auto Fx = F.pow(x);
    cout << print_to_string_i128(Fx.get_rlwe(0).get_poly(0).get_coeff(0)) << "\n";
    for (const auto& rlwe_el : Fx) {
        for (const auto& poly_el : rlwe_el) {
            for (const auto& coeff : poly_el)
                cout << print_to_string_i128(coeff) << ", ";
            cout << "\n";
        }
        cout << "\n\n\n\n";
    }
    return 0;
}