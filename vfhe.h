// TODO should default constructors init T arr to nullptrs? size=0...
// TODO replace vectori128 by an array1d i128 subclass which:
    // TODO implements scalar mult (should be covered by base class anyway...)
#pragma once

#include <iostream>
#include <cassert>
#include <random>
#include <chrono>
#include "omp.h"
#include "shared.h"
#include "atcoder/convolution.hpp" // Download from AtCoder ACL

using namespace atcoder;

#ifdef TIMING_ON
#  define TIMING(x) x
#else
#  define TIMING(x)
#endif

#ifdef DEBUG_ON
#  define DEBUG(x) x
#else
#  define DEBUG(x)
#endif

#ifdef DEBUG1_ON
#  define DEBUG1(x) x
#else
#  define DEBUG1(x)
#endif

#ifndef N_THREADS
#  define N_THREADS 8
#endif

typedef struct {
    int calls_convolve[N_THREADS] = { 0 };
    i128 iter_ = 0;
    std::chrono::duration<double, std::milli> elapsed_verify{};
    std::chrono::duration<double, std::milli> elapsed_proof{};
    std::chrono::duration<double, std::milli> elapsed_controller{};
    std::chrono::duration<double, std::milli> elapsed_plant{};
    std::chrono::duration<double, std::milli> elapsed_total{};
    std::chrono::duration<double, std::milli> convolve[N_THREADS]{};
} times_and_counts;

// NOTE inline keyword for structs allows the struct to be used in multiple
// translation units without causing linker errors
inline times_and_counts times_counts = { 0 };

vector_i128 sample_discrete_gaussian(size_t N, double mu = 3.2, double sigma = 19.2) {
    vector_i128 result(N);
    std::random_device rd;
    std::mt19937 gen(rd());
    std::normal_distribution<double> dist(mu, sigma);
    for (size_t i = 0; i < N; i++) {
        i128 res = static_cast<i128>(std::round(dist(gen)));
        result[i] = mod_(res, FIELD_MODULUS);
    }
    #ifdef DEBUG1_ON
        for (size_t i = 0; i < N; i++)
            result.at(i) = 1;
    #endif
    return result;
}
vector_i128 sample_secret_key(size_t N) {
    // Sample a secret key for the RGSW scheme.
    // Each entry is -1, 0, or 1, with probabilities 0.25, 0.5, 0.25 respectively.
    vector_i128 s(N);
    std::random_device rd;
    std::mt19937 gen(rd());
    std::discrete_distribution<int> dist({0.25, 0.5, 0.25});
    for (size_t i = 0; i < N; ++i) {
        int val = dist(gen);
        // TODO uncomment and deal with potential check_val_bounds errors
        // if (val == 0) s[i] = -1;
        if (val == 0) s[i] = 1;
        else if (val == 1) s[i] = 0;
        else s[i] = 1;
    }
    #ifdef DEBUG1_ON
    for (size_t i = 0; i < N; i++)
        s.at(i) = 1;
    #endif
    return s;
}
// Sample a random polynomial of degree N-1 with coefficients in the range [0, q).
vector_i128 sample_random_polynomial(size_t N, i128 q) {
    vector_i128 poly(N);
    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_int_distribution<i128> dist(0, q - 1);
    for (size_t i = 0; i < N; ++i) {
        poly[i] = dist(gen);
    }
    #ifdef DEBUG1_ON
        for (size_t i = 0; i < N; i++)
            poly.at(i) = 1;
    #endif
    return poly;
}

// Sample a noise polynomial of degree N-1 with coefficients from the discrete Gaussian distribution.
vector_i128 sample_noise_polynomial(size_t N, double mu = 3.2, double sigma = 19.2) {
    return sample_discrete_gaussian(N, mu, sigma);
}
template<typename T, typename Derived>
class array1d {
private:
    size_t size_;
    T* arr;
public:
    array1d() : size_(0), arr(nullptr) {}
    array1d(size_t size) : size_(size), arr(new T[size]()) {
        arr = new T[size]();
    }
    // Copy constructor
    array1d(const array1d& other) : size_(other.size_), arr(new T[other.size_]) {
        for (size_t i = 0; i < size_; ++i) arr[i] = other.arr[i];
    }
    // Move constructor
    array1d(array1d&& other) noexcept : size_(other.size_), arr(other.arr) {
        other.size_ = 0;
        other.arr = nullptr;
    }
    ~array1d() {
        delete[] arr;
    }
    void check_index_bounds(size_t n) const {
        if (n >= size_) {
            throw std::out_of_range(
                "Index error: accessing arr[" + std::to_string(n) + "]"
                + " in a " + std::to_string(size_) + " element array."
            );
        }
    }
    // T& get(size_t n, bool disable_value_check = true) const {
    T& get(size_t n) const {
        DEBUG(check_index_bounds(n);)
        // if (!disable_value_check)
        DEBUG(check_value_bounds(arr[n]);)
        return arr[n];
    }
    void check_value_bounds(const T& val) const {
        if constexpr (std::is_same_v<T, i128>) {
            i128 min_val = 0;
            i128 max_val = GROUP_MODULUS - 1;
            if (val < min_val || val > max_val) {
            // if (val > max_val) {
                throw std::out_of_range(
                    "(array1d) Value out of range: " + print_to_string_i128(val)
                );
            }
        }
    }
    // void set(int n, T val, bool disable_value_check = false) { // FIXME should be T& ?
    void set(int n, const T& val) {
        DEBUG(check_index_bounds(n);)
        // if (!disable_value_check)
        DEBUG(check_value_bounds(val);)
        arr[n] = val;
    }
    size_t size() const {
        return size_;
    }
    // Copy assignment operator
    array1d& operator=(const array1d& other) {
        if (this != &other) {
            delete[] arr;
            size_ = other.size_;
            arr = new T[size_];
            for (size_t i = 0; i < size_; ++i) arr[i] = other.arr[i];
        }
        return *this;
    }
    // Move assignment operator
    array1d& operator=(array1d&& other) noexcept {
        if (this != &other) {
            delete[] arr;
            size_ = other.size_;
            arr = other.arr;
            other.size_ = 0;
            other.arr = nullptr;
        }
        return *this;
    }
    Derived operator*(const Derived& other) const {
        size_t N = size();
        Derived res(N);
        for (size_t i = 0; i < N; i++) {
            auto val = (get(i) * other.get(i));
            if constexpr (std::is_same_v<T, i128>)
                val = mod_(val, FIELD_MODULUS);
            res.set(i, val);
        }
        return res;
    }
    Derived operator*(const i128& scalar) const {
        size_t N = size();
        Derived res(N);
        for (size_t i = 0; i < N; i++) {
            auto val = get(i) * scalar;
            if constexpr (std::is_same_v<T, i128>)
                val = mod_(val, FIELD_MODULUS);
            res.set(i, val);
        }
        return res;
    }
    i128 group_mult_(i128 a, i128 b) const {
        return mod_(a * b, GROUP_MODULUS);
    }
    Derived group_mult(const Derived& other) const {
        size_t N = size();
        Derived res(N);
        for (size_t i = 0; i < N; i++) {
            T val;
            if constexpr (std::is_same_v<T, i128>) {
                // val = (get(i) * other.get(i));
                // val = mod_(val, GROUP_MODULUS);
                val = group_mult_(get(i), other.get(i));
            } else {
                val = get(i).group_mult(other.get(i));
            }
            res.set(i, val);
        }
        return res;
    }
    Derived group_mult(const i128& scalar) const {
        size_t N = size();
        Derived res(N);
        for (size_t i = 0; i < N; i++) {
            T val;
            if constexpr (std::is_same_v<T, i128>) {
                // val = get(i) * scalar;
                // val = mod_(val, GROUP_MODULUS);
                val = group_mult_(get(i) * scalar, GROUP_MODULUS);
            } else {
                val = get(i).group_mult(scalar);
            }
            res.set(i, val);
        }
        return res;
    }
    Derived operator+(const Derived& other) const {
        size_t N = size();
        Derived res(N);
        for (size_t i = 0; i < N; i++) {
            auto val = get(i) + other.get(i);
            if constexpr (std::is_same_v<T, i128>)
                val = mod_(val, FIELD_MODULUS);
            res.set(i, val);
        }
        return res;
    }

    Derived operator-(const Derived& other) const {
        size_t N = size();
        Derived neg_other(N);
        for (size_t i = 0; i < N; i++) {
            auto val = get(i) - other.get(i);
            if constexpr (std::is_same_v<T, i128>)
                val = mod_(val, FIELD_MODULUS);
            neg_other.set(i, val);
        }
        return neg_other;
    }
    // raise self to other
    Derived pow(const Derived& other) const {
        size_t N = size();
        Derived res(N);
        if constexpr (std::is_same_v<T, i128>) {
            for (size_t i = 0; i < N; i++) {
                auto val = pow_(get(i), other.get(i));
                res.set(i, val);
            }
            return res;
        } else {
            for (size_t i = 0; i < N; i++) {
                auto val = get(i).pow(other.get(i));
                res.set(i, val);
            }
            return res;
        }
    }
    // raise scalar to self
    Derived pow(const i128 scalar) const {
        size_t N = size();
        Derived res(N);
        if constexpr (std::is_same_v<T, i128>) {
            for (size_t i = 0; i < N; i++) {
                auto val = pow_(scalar, get(i));
                res.set(i, val);
            }
            return res;
        } else {
            for (size_t i = 0; i < N; i++) {
                auto val = get(i).pow(scalar);
                res.set(i, val);
            }
            return res;
        }
    }
    // raise generator to self
    Derived pow() const {
        size_t N = size();
        Derived res(N);
        if constexpr (std::is_same_v<T, i128>) {
            for (size_t i = 0; i < N; i++) {
                auto val = pow_(GENERATOR, get(i));
                res.set(i, val);
            }
            return res;
        } else {
            for (size_t i = 0; i < N; i++) {
                auto val = get(i).pow();
                res.set(i, val);
            }
            return res;
        }
    }
    // base case binary modular exponentiation
    i128 pow_(i128 base, i128 power) const {
        power = mod_(power, FIELD_MODULUS);
        base = mod_(base, GROUP_MODULUS);
        i128 result = 1;
        while (power > 0) {
            bool is_power_odd = (power % 2) == 1;
            if (is_power_odd)
                result = (result * base) % GROUP_MODULUS;
            power >>= 1;
            base = (base * base) % GROUP_MODULUS;
        }
        return result;
    }

    T* begin() { return arr; }
    T* end() { return arr + size_; }
    const T* begin() const { return arr; }
    const T* end() const { return arr + size_; }

    void print() const {
        if constexpr (std::is_same_v<T, i128>) {
            std::cout << "{";
            for (size_t i = 0; i < size() - 1; i++) {
                std::cout << print_to_string_i128(get(i));
                std::cout << ", ";
            }
            std::cout << print_to_string_i128(get(size() - 1)) << "}";
        } else {
            std::cout << "{";
            for (size_t i = 0; i < size() - 1; i++) {
                get(i).print();
                std::cout << ", ";
            }
            get(size() - 1).print();
            std::cout << "}  ";
            std::cout << "\n";
        }
    }
};

template<typename T>
class array2d {
private:
    size_t rows_;
    size_t cols_;
    T* data;
    T** arr;
public:
    array2d() : rows_(0), cols_(0), data(nullptr), arr(nullptr) {}
    array2d(size_t rows, size_t cols) : rows_(rows), cols_(cols) {
        data = new T[rows * cols]();
        arr = new T*[rows];
        for (size_t i = 0; i < rows; i++) {
            arr[i] = &data[i * cols];
        }
    }
    ~array2d() {}
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
        DEBUG(check_index_bounds(row, col);)
        DEBUG(check_value_bounds(arr[row][col]);)
        return arr[row][col];
    }

    void check_value_bounds(const T& val) const {
        if constexpr (std::is_same_v<T, i128>) {
            i128 min_val = -1;
            min_val <<= 63;
            i128 max_val = (1UL << 63) - 1;
            if (val < min_val || val > max_val) {
                throw std::out_of_range(
                    "(array2d) Value out of range: " + print_to_string_i128(val)
                );
            }
        }
    }
    void set(size_t row, size_t col, const T& val) {
        DEBUG(check_index_bounds(row, col);)
        DEBUG(check_value_bounds(val);)
        arr[row][col] = val;
    }
    size_t n_rows() const {
        return rows_;
    }
    size_t n_cols() const {
        return cols_;
    }
    void print_i128() const {
        for (size_t i = 0; i < rows_; i++) {
            std::cout << "{";
            for (size_t j = 0; j < cols_ - 1; j++) {
                std::cout << print_to_string_i128(arr[i][j]) << ", ";
            }
            std::cout << print_to_string_i128(arr[i][cols_ - 1]) << "}\n";
        }
    }
    void print_array1d() const {
        for (size_t i = 0; i < rows_; i++) {
            for (size_t j = 0; j < cols_ - 1; j++) {
                std::cout << "[" << i << "]" << "[" << j << "]:\n";
                arr[i][j].print();
            }
        }
    }
    void print() const {
        if constexpr (std::is_same_v<T, i128>) {
            print_i128();
        } else {
            print_array1d();
        }
    }
};

class poly : public array1d<i128, poly> {
private:
    bool isNTT = false;
public:
    poly() : array1d<i128, poly>() {}
    poly(size_t N) : array1d<i128, poly>(N) {
        for (size_t i = 0; i < N; i++)
            set(i, 0);
    }
    auto& get_coeff(size_t n) const {
        return get(n);
    }
    size_t n_coeffs() const {
        return size();
    }
    auto get_hash(vector_i128 eval_pows) const {
        i128 hash = 0;
        for (size_t i = 0; i < size(); i++) {
            hash = mod_(hash + get(i) * eval_pows.at(i), FIELD_MODULUS);
        }
        return hash;
    }
    // See inline definition below class hashed_a_poly (avoids circular definitions)
    auto get_hash_a(vector_i128 eval_pows) const;

    std::vector<i128> convolution_(const std::vector<i128>& a, const std::vector<i128>& b) const {
        i128 n = a.size();
        std::vector<i128> a_pad = a, b_pad = b;
        std::vector<i128> conv = convolution<FIELD_MODULUS>(a_pad, b_pad);
        assert(conv.size() == 2 * n - 1);
        return conv;
    }

    std::vector<i128> conv_to_nega_(const std::vector<i128> conv) const {
        i128 n = (conv.size() + 1) / 2;
        std::vector<i128> res(n);
        for (i128 i = 0; i < n - 1; i++) {
            i128 val = mod_(conv[i] - conv[i + n], FIELD_MODULUS);
            res[i] = val;
        }
        res.at(n - 1) = conv.at(n - 1);
        return res;
    }
    std::vector<i128> negacyclic_convolution_(const std::vector<i128>& a, const std::vector<i128>& b) const {
        std::vector<i128> res = convolution_(a, b);
        std::vector<i128> res1 = conv_to_nega_(res);
        return res1;
    }
    poly convolve_naive(const poly& other) const {
        poly convolved(2 * n_coeffs() - 1);
        for (size_t i = 0; i < n_coeffs(); i++) {
            for (size_t j = 0; j < n_coeffs(); j++) {
                i128 val = mod_(get_coeff(i) * other.get_coeff(j), FIELD_MODULUS);
                val = mod_(convolved.get(i + j) + val, FIELD_MODULUS);
                convolved.set(i + j, val);
            }
        }
        return convolved;
    }
    poly convolve_ntt(const poly& other) const {
        TIMING(auto start = std::chrono::high_resolution_clock::now();)
        // TIMING(int thread_num = omp_get_thread_num();)
        TIMING(int thread_num = 0;)
        TIMING(times_counts.calls_convolve[thread_num] += 1;)

        // turn *this and other into vector_i128's
        size_t n = n_coeffs();
        vector_i128 a(n);
        vector_i128 b(n);
        for (size_t i = 0; i < n; i++) {
            a.at(i) = get(i);
            b.at(i) = other.get(i);
        }
        // call negacyclic_convolution_()
        vector_i128 conv = convolution_(a, b);
        size_t n_conv = 2 * n - 1;
        assert(conv.size() == n_conv);
        // create poly of result
        poly res(n_conv);
        for (size_t i = 0; i < n_conv; i++)
            res.set(i, conv.at(i));

        TIMING(auto end = std::chrono::high_resolution_clock::now();)
        TIMING(times_counts.convolve[thread_num] += end - start;)
        return res;
    }
    poly convolve(const poly& other) const {
        assert(n_coeffs() == other.n_coeffs());
        bool using_ntt = true;
        if (using_ntt)
            return convolve_ntt(other);
        else
            return convolve_naive(other);
    }
    poly nega_ntt(const poly& other) const {
        // turn *this and other into vector_i128's
        size_t n = n_coeffs();
        vector_i128 a(n);
        vector_i128 b(n);
        for (size_t i = 0; i < n; i++) {
            a.at(i) = get(i);
            b.at(i) = other.get(i);
        }
        // call negacyclic_convolution_()
        vector_i128 conv = negacyclic_convolution_(a, b);
        assert(conv.size() == n);
        // create poly of result
        poly res(n);
        for (size_t i = 0; i < n; i++)
            res.set(i, conv.at(i));
        return res;
    }
    poly nega_naive(const poly& other) const {
        size_t N = n_coeffs();
        poly neg_conv(N);
        for (size_t i = 0; i < N; i++) {
            for (size_t j = 0; j < i + 1; j++) {
                i128 val = mod_(get_coeff(j) * other.get_coeff(i - j), FIELD_MODULUS);
                val = mod_(neg_conv.get(i) + val, FIELD_MODULUS);
                neg_conv.set(i, val);
            }
            for (size_t j = i + 1; j < N; j++) {
                i128 val = mod_(get_coeff(j) * other.get_coeff(N + i - j), FIELD_MODULUS);
                val = mod_(neg_conv.get(i) - val, FIELD_MODULUS);
                neg_conv.set(i, val);
            }
        }
        return neg_conv;
    }
    using array1d<i128, poly>::operator*;
    poly operator*(const poly& other) const {
        size_t N = n_coeffs();
        assert(N == other.n_coeffs());
        bool using_ntt = true;
        if (using_ntt)
            return nega_ntt(other);
        else
            return nega_naive(other);
    }
    auto conv_to_nega(size_t N) const {
        assert(n_coeffs() == 2 * N - 1);
        i128 conv_degree = n_coeffs() - 1;
        // HACK can't return *this after defining array1d move semantics
        // assert(conv_degree >= N);
        if (conv_degree < N)
            return *this;
        poly negacyclic(N);
        for (size_t i = 0; i < N; i++)
            negacyclic.set(i, get_coeff(i));
        for (size_t i = N; i < conv_degree + 1; i++) {
            i128 coeff = negacyclic.get_coeff(i - N) - get_coeff(i);
            negacyclic.set(i - N, mod_(coeff, FIELD_MODULUS));
        }
        return negacyclic;
    }
};

class hashed_a_poly : public array1d<i128, hashed_a_poly> {
private:
public:
    hashed_a_poly() : array1d<i128, hashed_a_poly>() {}
    hashed_a_poly(size_t n_hashed_a_coeffs) : array1d<i128, hashed_a_poly>(n_hashed_a_coeffs) {
        for (size_t i = 0; i < n_hashed_a_coeffs; i++)
            set(i, 0);
    }
    auto& get_hashed_a_coeff(size_t n) const {
        return get(n);
    }
    size_t n_hashed_a_coeffs() const {
        return size();
    }
    i128 get_hash_sec(const poly& other) const {
        // HACK get around gr having hashed_a_polys of length 2n-1, but x_nega_ only length n
        assert(size() == other.size() || size() == (2 * other.size() - 1));
        i128 result = 1;
        poly raised(other.size());
        for (size_t i = 0; i < other.size(); i++) {
            raised.set(i, pow_(get(i), other.get(i)));
        }
        for (size_t i = 0; i < other.size(); i++) {
            result = group_mult_(result, raised.get(i));
        }
        return result;
    }
};

inline auto poly::get_hash_a(vector_i128 eval_pows) const {
    auto hash = get_hash(eval_pows);
    auto N = size();
    hashed_a_poly ha_poly(N);
    // TODO optimise
    vector_i128 hashed_a_vec = scalar_vec_mult(hash, eval_pows, FIELD_MODULUS);
    for (size_t i = 0; i < N; i++) {
        auto val = hashed_a_vec.at(i);
        ha_poly.set(i, val);
    }
    return ha_poly;
}

class hashed_rlwe_decomp : public array1d<i128, hashed_rlwe_decomp> {
private:
public:
    // NOTE hashed_rlwe always has two hashed_polys
    hashed_rlwe_decomp() : array1d<i128, hashed_rlwe_decomp>() {}
    hashed_rlwe_decomp(size_t n_hashed_polys) : array1d<i128, hashed_rlwe_decomp>(n_hashed_polys) {}
    auto& get_hashed_poly(size_t n) const {
        return get(n);
    }
    size_t n_hashed_polys() const {
        return size();
    }
};

class rlwe_decomp : public array1d<poly, rlwe_decomp> {
private:
public:
    rlwe_decomp() : array1d<poly, rlwe_decomp>() {}
    rlwe_decomp(size_t n_polys) : array1d<poly, rlwe_decomp>(n_polys) {}
    rlwe_decomp(
        size_t n_polys, size_t n_coeffs
    ) : array1d<poly, rlwe_decomp>(n_polys) {
        for (size_t i = 0; i < n_polys; i++) {
            set(i, poly(n_coeffs));
        }
    }
    ~rlwe_decomp() {}
    poly& get_poly(size_t n) const {
        return get(n);
    }
    // calls poly's get_hash and returns hashed_rlwe_decomp
    auto get_hash(vector_i128 eval_pows) const {
        hashed_rlwe_decomp hash(size()); // FIXME n_polys, do better
        for (size_t i = 0; i < size(); i++) {
            hash.set(i, get(i).get_hash(eval_pows));
        }
        return hash;
    }
    size_t n_polys() const {
        return size();
    }
    size_t n_coeffs() const {
        return get_poly(0).size();
    }
};

class hashed_rlwe_decomp_vec : public array1d<hashed_rlwe_decomp, hashed_rlwe_decomp_vec> {
private:
public:
    hashed_rlwe_decomp_vec() : array1d<hashed_rlwe_decomp, hashed_rlwe_decomp_vec>() {}
    hashed_rlwe_decomp_vec(size_t n_hashed_rlwe_decomps) : array1d<hashed_rlwe_decomp, hashed_rlwe_decomp_vec>(n_hashed_rlwe_decomps) {}
    hashed_rlwe_decomp_vec(
        size_t n_hashed_rlwe_decomps, size_t n_hashed_polys
    ) : array1d<hashed_rlwe_decomp, hashed_rlwe_decomp_vec>(n_hashed_rlwe_decomps) {
        for (size_t i = 0; i < n_hashed_rlwe_decomps; i++)
            set(i, hashed_rlwe_decomp(n_hashed_polys));
    }
    auto& get_hashed_rlwe_decomp(size_t n) const {
        return get(n);
    }
    size_t n_hashed_rlwe_decomps() const {
        return size();
    }
    size_t n_hashed_polys() const {
        return get_hashed_rlwe_decomp(0).n_hashed_polys();
    }
};

class rlwe_decomp_vec : public array1d<rlwe_decomp, rlwe_decomp_vec> {
private:
public:
    rlwe_decomp_vec() : array1d<rlwe_decomp, rlwe_decomp_vec>() {}
    rlwe_decomp_vec(
        size_t n_rlwe_decomps
    ) : array1d<rlwe_decomp, rlwe_decomp_vec>(n_rlwe_decomps) {}
    rlwe_decomp_vec(
        size_t n_rlwe_decomps, size_t n_polys, size_t n_coeffs
    ) : array1d<rlwe_decomp, rlwe_decomp_vec>(n_rlwe_decomps) {
        for (size_t i = 0; i < n_rlwe_decomps; i++)
            // FIXME n_polys should match n_rlwes for rgsw. rethink terminology
            set(i, rlwe_decomp(n_polys, n_coeffs));
    }
    ~rlwe_decomp_vec() {}
    rlwe_decomp& get_rlwe_decomp(size_t n) const {
        return get(n);
    }
    size_t n_rlwe_decomps() const {
        return size();
    }
    size_t n_polys() const {
        return get_rlwe_decomp(0).size();
    }
    size_t n_coeffs() const {
        return get_rlwe_decomp(0).get_poly(0).size();
    }
    // calls rlwe_decomp's get_hash and returns hashed_rlwe_decomp_vec
    auto get_hash(vector_i128 eval_pows) const {
        hashed_rlwe_decomp_vec hash(n_rlwe_decomps(), get(0).size()); // FIXME n_polys, do better (global?)
        for (size_t i = 0; i < size(); i++) {
            hash.set(i, get(i).get_hash(eval_pows));
        }
        return hash;
    }
};

class hashed_a_rlwe : public array1d<hashed_a_poly, hashed_a_rlwe> {
private:
public:
    // NOTE hashed_a_rlwe always has two hashed_a_polys
    hashed_a_rlwe() : array1d<hashed_a_poly, hashed_a_rlwe>() {}
    hashed_a_rlwe(size_t n_hashed_a_polys) : array1d<hashed_a_poly, hashed_a_rlwe>(n_hashed_a_polys) {}
    hashed_a_rlwe(size_t n_hashed_a_polys, size_t n_coeffs) : array1d<hashed_a_poly, hashed_a_rlwe>(n_hashed_a_polys) {
        assert(n_hashed_a_polys == N_POLYS_IN_RLWE);
        for (size_t i = 0; i < n_hashed_a_polys; i++) {
            set(i, hashed_a_poly(n_coeffs));
        }
    }
    auto& get_hashed_a_poly(size_t n) const {
        return get(n);
    }
    size_t n_hashed_a_polys() const {
        return size();
    }
    // void set_coeffs_to_one() {
    //     for (size_t i = 0; i < size(); i++) {
    //         set(i, 1);
    //     }
    // }
};

class hashed_rlwe : public array1d<i128, hashed_rlwe> {
private:
public:
    // NOTE hashed_rlwe always has two hashed_polys
    hashed_rlwe() : array1d<i128, hashed_rlwe>() {}
    hashed_rlwe(size_t n_hashed_polys) : array1d<i128, hashed_rlwe>(n_hashed_polys) {
        assert(n_hashed_polys == N_POLYS_IN_RLWE);
    }
    auto& get_hashed_poly(size_t n) const {
        return get(n);
    }
    size_t n_hashed_polys() const {
        return size();
    }
    void set_coeffs_to_one() {
        for (size_t i = 0; i < size(); i++) {
            set(i, 1);
        }
    }
};

class rlwe : public array1d<poly, rlwe> {
private:
public:
    // NOTE rlwe always has two polys
    rlwe() : array1d<poly, rlwe>() {}
    rlwe(size_t n_polys) : array1d<poly, rlwe>(n_polys) {
        assert(n_polys == 2);
    }
    rlwe(size_t n_polys, size_t n_coeffs) : array1d<poly, rlwe>(n_polys) {
        assert(n_polys == 2);
        for (size_t i = 0; i < n_polys; i++) {
            set(i, poly(n_coeffs));
        }
    }
    auto& get_poly(size_t n) const {
        return get(n);
    }
    // TODO move all get_hash() to array1d?
    // calls poly's get_hash and returns hashed_rlwe
    auto get_hash(vector_i128 eval_pows) const {
        hashed_rlwe hash(N_POLYS_IN_RLWE); // FIXME n_polys, do better (global?)
        for (size_t i = 0; i < size(); i++) {
            hash.set(i, get(i).get_hash(eval_pows));
        }
        return hash;
    }
    // calls poly's get_hash and returns hashed_rlwe
    auto get_hash_a(vector_i128 eval_pows) const {
        hashed_a_rlwe hash_a(N_POLYS_IN_RLWE, eval_pows.size() / 2); // FIXME n_polys, do better (global?)
        for (size_t i = 0; i < size(); i++) {
            hash_a.set(i, get(i).get_hash_a(eval_pows));
        }
        return hash_a;
    }
    void set_coeffs_to_one() {
        for (auto& p : *this) {
            for (size_t j = 0; j < p.size(); j++) {
                p.set(j, 1);
            }
        }
    }
    size_t n_polys() const {
        return size();
    }
    size_t n_coeffs() const {
        return get_poly(0).size();
    }
    // creates d polynomials for each of the two polynomials in the rlwe object,
    // where the i'th coefficient has been decomposed with base v and depth d.
    // The i'th coefficient is spread across the i'th coefficients of the d polynomials.
    rlwe_decomp decompose(const i128& v_, const i128& d) const {
        // v - power of 2, s.t. v^{d-1} < q < v^d
        i128 power = static_cast<i128>(std::ceil((std::log2(FIELD_MODULUS) / d)));
        // TODO remove double up of computation of v in run_control_loop
        i128 v = 1 << power;
        assert(v == v_);
        assert(d >= 1);
        // FIXME not sure if we really need the lower bounds check. Fails sometimes
        // assert(v**(d-1) < q and q <= v**d)
        assert(FIELD_MODULUS <= static_cast<i128>(std::pow(v, d)));

        // decompose poly into d polynomials of degree d-1
        // polys = []
        auto pow_ = [](i128 base, i128 exp) {
            i128 res = 1;
            for (i128 i = 0; i < exp; ++i) {
                res *= base;
            }
            return res;
        };
        rlwe_decomp polys(2 * d, n_coeffs());
        for (size_t k = 0; k < N_POLYS_IN_RLWE; k++) {
            poly pol = get_poly(k);
            for (size_t i = 0; i < n_coeffs(); i++) {
                for (size_t j = 0; j < d; j++) {
                    i128 decomped_coeff = (v - 1) & (pol.get_coeff(i) / pow_(v, j));
                    assert(decomped_coeff < v);
                    polys.get_poly(d * k + j).set(i, decomped_coeff);
                }
            }
        }
        return polys;
    }
    auto conv_to_nega(size_t N) const {
        rlwe conv(N_POLYS_IN_RLWE);
        for (size_t i = 0; i < N_POLYS_IN_RLWE; i++) {
            poly p = get_poly(i).conv_to_nega(N);
            conv.set(i, p);
        }
        return conv;
    }
};

class hashed_rlwe_vec : public array1d<hashed_rlwe, hashed_rlwe_vec> {
private:
public:
    hashed_rlwe_vec() : array1d<hashed_rlwe, hashed_rlwe_vec>() {}
    hashed_rlwe_vec(size_t n_hashed_rlwes) : array1d<hashed_rlwe, hashed_rlwe_vec>(n_hashed_rlwes) {}
    hashed_rlwe_vec(
        size_t n_hashed_rlwes, size_t n_hashed_polys
    ) : array1d<hashed_rlwe, hashed_rlwe_vec>(n_hashed_rlwes) {
        for (size_t i = 0; i < n_hashed_rlwes; i++) {
            assert(n_hashed_polys == 2);
            set(i, hashed_rlwe(n_hashed_polys));
        }
    }
    auto& get_hashed_rlwe(size_t n) const {
        return get(n);
    }
    size_t n_hashed_rlwes() const {
        return size();
    }
    size_t n_hashed_polys() const {
        return get_hashed_rlwe(0).n_hashed_polys();
    }
};

class rlwe_vec : public array1d<rlwe, rlwe_vec> {
private:
public:
    rlwe_vec() : array1d<rlwe, rlwe_vec>() {}
    rlwe_vec(size_t n_rlwes) : array1d<rlwe, rlwe_vec>(n_rlwes) {}
    rlwe_vec(
        size_t n_rlwes, size_t n_polys, size_t n_coeffs
    ) : array1d<rlwe, rlwe_vec>(n_rlwes) {
        for (size_t i = 0; i < n_rlwes; i++)
            set(i, rlwe(n_polys, n_coeffs));
    }
    auto& get_rlwe(size_t n) const {
        return get(n);
    }
    // TODO move all get_hash() to array1d?
    // calls rlwe's get_hash and returns hashed_rlwe_vec
    auto get_hash(vector_i128 eval_pows) const {
        hashed_rlwe_vec hash(size(), get(0).size()); // FIXME n_polys, do better (global?)
        for (size_t i = 0; i < size(); i++) {
            hash.set(i, get(i).get_hash(eval_pows));
        }
        return hash;
    }
    size_t n_rlwes() const {
        return size();
    }
    size_t n_polys() const {
        return get_rlwe(0).size();
    }
    size_t n_coeffs() const {
        return get_rlwe(0).get_poly(0).size();
    }
    // calls rlwe's decompose on each rlwe element and returns a rlwe_decomp_vec
    auto decompose(const i128& v, const i128& d) const {
        rlwe_decomp_vec decomps(size(), 2 * d, n_coeffs());
        for (size_t i = 0; i < size(); i++)
            decomps.set(i, get_rlwe(i).decompose(v, d));
        return decomps;
    }
    auto conv_to_nega(size_t N) const {
        rlwe_vec conv(n_rlwes());
        for (size_t i = 0; i < n_rlwes(); i++) {
            rlwe r = get_rlwe(i).conv_to_nega(N);
            conv.set(i, r);
        }
        return conv;
    }
};

class hashed_a_rgsw : public array1d<hashed_a_rlwe, hashed_a_rgsw> {
private:
public:
    hashed_a_rgsw() : array1d<hashed_a_rlwe, hashed_a_rgsw>() {}
    hashed_a_rgsw(size_t n_hashed_a_rlwes) : array1d<hashed_a_rlwe, hashed_a_rgsw>(n_hashed_a_rlwes) {
    }
    hashed_a_rgsw(
        size_t n_hashed_a_rlwes, size_t n_hashed_a_polys, size_t n_coeffs
    ) : array1d<hashed_a_rlwe, hashed_a_rgsw>(n_hashed_a_rlwes) {
        for (size_t i = 0; i < n_hashed_a_rlwes; i++)
            set(i, hashed_a_rlwe(n_hashed_a_polys, n_coeffs));
    }
    auto& get_hashed_a_rlwe(size_t n) const {
        return get(n);
    }
    size_t n_hashed_a_rlwes() const {
        return size();
    }
    size_t n_hashed_a_polys() const {
        return get_hashed_a_rlwe(0).n_hashed_a_polys();
    }
    // using array1d<hashed_a_rlwe, hashed_a_rgsw>::pow;
    hashed_rlwe get_hash_sec(const rlwe_decomp& other) const {
        size_t n_polys_ = n_hashed_a_polys();
        assert(n_polys_ == N_POLYS_IN_RLWE);
        size_t N = size();
        assert(N == other.size());
        hashed_rlwe_vec res_vec(N, n_polys_);
        // #pragma omp parallel for num_threads(2)
        for (size_t i = 0; i < N; i++) {
            auto val0 = get(i).get(0).get_hash_sec(other.get(i));
            auto val1 = get(i).get(1).get_hash_sec(other.get(i));
            hashed_rlwe& res = res_vec.get(i);
            res.set(0, val0);
            res.set(1, val1);
        }

        hashed_rlwe res(n_polys_);
        res.set_coeffs_to_one();
        for (auto& hashed_a_rlwe_ : res_vec) {
            res = res.group_mult(hashed_a_rlwe_);
        }
        return res;
    }
};

class hashed_rgsw : public array1d<hashed_rlwe, hashed_rgsw> {
private:
public:
    hashed_rgsw() : array1d<hashed_rlwe, hashed_rgsw>() {}
    hashed_rgsw(size_t n_hashed_rlwes) : array1d<hashed_rlwe, hashed_rgsw>(n_hashed_rlwes) {
    }
    hashed_rgsw(
        size_t n_hashed_rlwes, size_t n_hashed_polys
    ) : array1d<hashed_rlwe, hashed_rgsw>(n_hashed_rlwes) {
        for (size_t i = 0; i < n_hashed_rlwes; i++)
            set(i, hashed_rlwe(n_hashed_polys));
    }
    auto& get_hashed_rlwe(size_t n) const {
        return get(n);
    }
    size_t n_hashed_rlwes() const {
        return size();
    }
    size_t n_hashed_polys() const {
        return get_hashed_rlwe(0).n_hashed_polys();
    }
    using array1d<hashed_rlwe, hashed_rgsw>::pow;
    hashed_rlwe pow(const hashed_rlwe_decomp& other) const {
        size_t n_polys_ = n_hashed_polys();
        assert(n_polys_ == N_POLYS_IN_RLWE);
        size_t N = size();
        assert(N == other.size());
        hashed_rlwe_vec res_vec(N, n_polys_);
        // #pragma omp parallel for num_threads(2)
        for (size_t i = 0; i < N; i++) {
            // FIXME a bit hacky - obj.pow(scalar) is scalar^obj, do we need obj^scalar???
            auto val0 = pow_(get(i).get(0), (other.get(i)));
            auto val1 = pow_(get(i).get(1), (other.get(i)));
            hashed_rlwe& res = res_vec.get(i);
            res.set(0, val0);
            res.set(1, val1);
        }

        hashed_rlwe res(n_polys_);
        res.set_coeffs_to_one();
        for (auto& hashed_rlwe_ : res_vec) {
            res = res.group_mult(hashed_rlwe_);
        }
        return res;
    }
};

class rgsw : public array1d<rlwe, rgsw> {
private:
public:
    rgsw() : array1d<rlwe, rgsw>() {}
    rgsw(size_t n_rlwes) : array1d<rlwe, rgsw>(n_rlwes) {}
    rgsw(
        size_t n_rlwes, size_t n_polys, size_t n_coeffs
    ) : array1d<rlwe, rgsw>(n_rlwes) {
        for (size_t i = 0; i < n_rlwes; i++)
            set(i, rlwe(n_polys, n_coeffs));
    }
    ~rgsw() {}
    rlwe& get_rlwe(size_t n) const {
        return get(n);
    }
    size_t n_rlwes() const {
        return size();
    }
    size_t n_polys() const {
        return get_rlwe(0).n_polys();
    }
    size_t n_coeffs() const {
        return get_rlwe(0).get_poly(0).n_coeffs();
    }
    using array1d<rlwe, rgsw>::operator*;
    rlwe operator*(const rlwe_decomp& other) const {
        size_t N = size();
        assert(N == other.size());
        rlwe res(n_polys(), n_coeffs());
        // TODO wrap in loop
        poly& p0 = res.get(0);
        poly& p1 = res.get(1);
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

    rlwe convolve(const rlwe_decomp& other) const {
        size_t N = size();
        assert(N == other.size());
        rlwe res(n_polys(), 2 * n_coeffs() - 1);
        // TODO wrap in loop
        poly& p0 = res.get(0);
        poly& p1 = res.get(1);
        // TODO wrap in outer loop
        // FIXME need thread array to accumulate sum
        // #pragma omp parallel for num_threads(N_THREADS)
        for (size_t i = 0; i < N; i++) {
            auto val0 = get(i).get(0).convolve(other.get(i));
            auto val1 = get(i).get(1).convolve(other.get(i));
            p0 = p0 + val0;
            p1 = p1 + val1;
        }
        res.set(0, p0);
        res.set(1, p1);
        return res;
    }

    using array1d<rlwe, rgsw>::pow;
    rlwe pow(const rlwe_decomp& other) const {
        size_t n_polys_ = n_polys();
        size_t n_coeffs_ = n_coeffs();
        size_t N = size();
        assert(N == other.size());
        rlwe_vec res_vec(N, n_polys_, n_coeffs_);
        // TODO wrap in loop
        // #pragma omp parallel for num_threads(N_THREADS)
        for (size_t i = 0; i < N; i++) {
            rlwe& res = res_vec.get(i);
            auto val0 = get(i).get(0).pow(other.get(i));
            auto val1 = get(i).get(1).pow(other.get(i));
            res.set(0, val0);
            res.set(1, val1);
        }

        rlwe res(n_polys_, n_coeffs_);
        res.set_coeffs_to_one();
        for (auto& rlwe_ : res_vec) {
            res = res.group_mult(rlwe_);
        }
        return res;
    }
    // calls rlwe's get_hash for all the rlwe's in (*this)
    auto get_hash(vector_i128 eval_pows) const {
        hashed_rgsw hash(size(), n_polys()); // FIXME n_polys, do better (global?)
        for (size_t i = 0; i < size(); i++) {
            hash.set(i, get(i).get_hash(eval_pows));
        }
        return hash;
    }
    // calls rlwe's get_hash_a for all the rlwe's in (*this)
    auto get_hash_a(vector_i128 eval_pows) const {
        hashed_a_rgsw hash(size(), n_polys(), n_coeffs()); // FIXME n_polys, do better (global?)
        for (size_t i = 0; i < size(); i++) {
            hash.set(i, get(i).get_hash_a(eval_pows));
        }
        return hash;
    }
};

class hashed_rgsw_vec : public array1d<hashed_rgsw, hashed_rgsw_vec> {
private:
public:
    hashed_rgsw_vec() : array1d<hashed_rgsw, hashed_rgsw_vec>() {}
    hashed_rgsw_vec(size_t n_hashed_rgsws) : array1d<hashed_rgsw, hashed_rgsw_vec>(n_hashed_rgsws) {}
    hashed_rgsw_vec(
        size_t n_hashed_rgsws, size_t n_hashed_rlwes, size_t n_hashed_polys
    ) : array1d<hashed_rgsw, hashed_rgsw_vec>(n_hashed_rgsws) {
        for (size_t i = 0; i < n_hashed_rgsws; i++)
            set(i, hashed_rgsw(n_hashed_rlwes, n_hashed_polys));
    }
    ~hashed_rgsw_vec() {}
    hashed_rgsw& get_hashed_rgsw(size_t n) const {
        return get(n);
    }
    size_t n_hashed_rgsws() const {
        return size();
    }
    size_t n_hashed_rlwes() const {
        size_t n_hashed_rlwes = get_hashed_rgsw(0).n_hashed_rlwes();
        return n_hashed_rlwes;
    }
    size_t n_hashed_polys() const {
        size_t n_hashed_polys = get_hashed_rgsw(0).n_hashed_polys();
        assert(n_hashed_polys == 2);
        return n_hashed_polys;
    }
    using array1d<hashed_rgsw, hashed_rgsw_vec>::pow;
    hashed_rlwe pow(const hashed_rlwe_decomp_vec& other) const {
        size_t n = other.size();
        assert(size() == n);
        size_t n_hashed_polys_ = n_hashed_polys();
        assert(n_hashed_polys_ == 2);
        hashed_rlwe_vec res_vec(n, n_hashed_polys_);
        hashed_rlwe res(n_hashed_polys_);
        res.set_coeffs_to_one();
        // #pragma omp parallel for num_threads(12)
        for (size_t i = 0; i < n; i++) {
            auto val = get(i).pow(other.get(i));
            res_vec.set(i, val);
        }
        for (auto& prod : res_vec) {
            res = res.group_mult(prod);
        }
        return res;
    }
};

class hashed_a_rgsw_vec : public array1d<hashed_a_rgsw, hashed_a_rgsw_vec> {
private:
public:
    hashed_a_rgsw_vec() : array1d<hashed_a_rgsw, hashed_a_rgsw_vec>() {}
    hashed_a_rgsw_vec(size_t n_hashed_a_rgsws) : array1d<hashed_a_rgsw, hashed_a_rgsw_vec>(n_hashed_a_rgsws) {}
    hashed_a_rgsw_vec(
        size_t n_hashed_a_rgsws, size_t n_hashed_a_rlwes, size_t n_hashed_a_polys, size_t n_coeffs
    ) : array1d<hashed_a_rgsw, hashed_a_rgsw_vec>(n_hashed_a_rgsws) {
        for (size_t i = 0; i < n_hashed_a_rgsws; i++)
            set(i, hashed_a_rgsw(n_hashed_a_rlwes, n_hashed_a_polys, n_coeffs));
    }
    ~hashed_a_rgsw_vec() {}
    hashed_a_rgsw& get_hashed_a_rgsw(size_t n) const {
        return get(n);
    }
    size_t n_hashed_a_rgsws() const {
        return size();
    }
    size_t n_hashed_a_rlwes() const {
        size_t n_hashed_a_rlwes = get_hashed_a_rgsw(0).n_hashed_a_rlwes();
        return n_hashed_a_rlwes;
    }
    size_t n_hashed_a_polys() const {
        size_t n_hashed_a_polys = get_hashed_a_rgsw(0).n_hashed_a_polys();
        assert(n_hashed_a_polys == 2);
        return n_hashed_a_polys;
    }
    // using array1d<hashed_a_rgsw, hashed_a_rgsw_vec>::pow;
    hashed_rlwe get_hash_sec(const rlwe_decomp_vec& other) const {
        size_t n = other.size();
        assert(size() == n);
        size_t n_hashed_a_polys_ = n_hashed_a_polys();
        assert(n_hashed_a_polys_ == 2);
        hashed_rlwe_vec res_vec(n, n_hashed_a_polys_);
        hashed_rlwe res(n_hashed_a_polys_);
        res.set_coeffs_to_one();
        // #pragma omp parallel for num_threads(12)
        for (size_t i = 0; i < n; i++) {
            auto val = get(i).get_hash_sec(other.get(i));
            res_vec.set(i, val);
        }
        for (auto& prod : res_vec) {
            res = res.group_mult(prod);
        }
        return res;
    }
};


class rgsw_vec : public array1d<rgsw, rgsw_vec> {
private:
public:
    rgsw_vec() : array1d<rgsw, rgsw_vec>() {}
    rgsw_vec(size_t n_rgsws) : array1d<rgsw, rgsw_vec>(n_rgsws) {}
    rgsw_vec(
        size_t n_rgsws, size_t n_rlwes, size_t n_polys, size_t n_coeffs
    ) : array1d<rgsw, rgsw_vec>(n_rgsws) {
        for (size_t i = 0; i < n_rgsws; i++)
            set(i, rgsw(n_rlwes, n_polys, n_coeffs));
    }
    ~rgsw_vec() {}
    rgsw& get_rgsw(size_t n) const {
        return get(n);
    }
    size_t n_rgsws() const {
        return size();
    }
    size_t n_rlwes() const {
        return get_rgsw(0).size();
    }
    size_t n_polys() const {
        return get_rgsw(0).n_polys();
    }
    size_t n_coeffs() const {
        return get_rgsw(0).n_coeffs();
    }
    using array1d<rgsw, rgsw_vec>::operator*;
    rlwe operator*(const rlwe_decomp_vec& other) const {
        size_t n = other.size();
        assert(size() == n);
        rlwe sum(n_polys(), n_coeffs());
        for (size_t i = 0; i < n; i++) {
            auto val = get(i) * other.get(i);
            sum = sum + val;
        }
        return sum;
    }
    using array1d<rgsw, rgsw_vec>::pow;
    rlwe pow(const rlwe_decomp_vec& other) const {
        size_t n = other.size();
        assert(size() == n);
        size_t n_polys_ = n_polys();
        size_t n_coeffs_ = n_coeffs();
        rlwe_vec res_vec(n, n_polys_, n_coeffs_);
        rlwe res(n_polys_, n_coeffs_);
        res.set_coeffs_to_one();
        // #pragma omp parallel for num_threads(N_THREADS)
        for (size_t i = 0; i < n; i++) {
            auto val = get(i).pow(other.get(i));
            res_vec.set(i, val);
        }
        for (auto& prod : res_vec) {
            res = res.group_mult(prod);
        }
        return res;
    }
    // calls rgsw's get_hash
    auto get_hash(vector_i128 eval_pows) const {
        hashed_rgsw_vec hash(size(), n_rlwes(), n_polys()); // FIXME n_polys, do better (global?)
        for (size_t i = 0; i < size(); i++) {
            hash.set(i, get(i).get_hash(eval_pows));
        }
        return hash;
    }
    // calls rgsw's get_hash_a
    auto get_hash_a(vector_i128 eval_pows) const {
        hashed_a_rgsw_vec hash(size(), n_rlwes(), n_polys(), n_coeffs()); // FIXME n_polys, do better (global?)
        for (size_t i = 0; i < size(); i++) {
            hash.set(i, get(i).get_hash_a(eval_pows));
        }
        return hash;
    }
};

class rgsw_mat : public array2d<rgsw> {
private:
public:
    rgsw_mat() : array2d<rgsw>() {}
    rgsw_mat(size_t rows, size_t cols) : array2d<rgsw>(rows, cols) {}
    rgsw_mat(
        size_t rows, size_t cols, size_t n_rlwes, size_t n_polys,
        size_t n_coeffs
    ) : array2d<rgsw>(rows, cols) {
        for (size_t i = 0; i < rows; i++) {
            for (size_t j = 0; j < cols; j++) {
                set(i, j, rgsw(n_rlwes, n_polys, n_coeffs));
            }
        }
    }

    rlwe_vec operator*(const rlwe_decomp_vec& other) const {
        size_t rows = n_rows();
        size_t cols = n_cols();
        assert(cols == other.size());
        size_t n_polys_ = n_polys();
        size_t n_coeffs_ = n_coeffs();
        size_t n_coeffs_other = other.n_coeffs();
        assert(n_coeffs_ == n_coeffs_other);

        rlwe_vec res(rows, n_polys_, n_coeffs_);
        for (size_t i = 0; i < rows; i++) {
            rlwe sum(n_polys_, n_coeffs_);
            for (size_t j = 0; j < cols; j++) {
                auto val = get(i, j) * other.get(j);
                sum = sum + val;
            }
            res.set(i, sum);
        }
        return res;
    }
    rlwe_vec convolve(const rlwe_decomp_vec& other) const {
        size_t rows = n_rows();
        size_t cols = n_cols();
        assert(cols == other.size());
        size_t n_polys_ = n_polys();
        size_t n_coeffs_ = n_coeffs();
        size_t n_coeffs_other = other.n_coeffs();
        assert(n_coeffs_ == n_coeffs_other);

        rlwe_vec res(rows, n_polys_, 2 * n_coeffs_ - 1);
        for (size_t i = 0; i < rows; i++) {
            rlwe sum(n_polys_, 2 * n_coeffs_ - 1);
            for (size_t j = 0; j < cols; j++) {
                auto val = get(i, j).convolve(other.get(j));
                sum = sum + val;
            }
            res.set(i, sum);
        }
        return res;
    }
    rlwe_vec pow(const rlwe_decomp_vec& other) const {
        size_t rows = n_rows();
        size_t cols = n_cols();
        assert(cols == other.size());
        rgsw& rg = get(0, 0);
        size_t n_coeffs = rg.n_coeffs();
        size_t n_polys = rg.n_polys();
        rlwe_vec res(rows, n_polys, n_coeffs);
        for (size_t i = 0; i < rows; i++) {
            rlwe sum(n_polys, n_coeffs);
            sum.set_coeffs_to_one();
            for (size_t j = 0; j < cols; j++) {
                auto val = get(i, j).pow(other.get(j));
                sum = sum.group_mult(val);
            }
            res.set(i, sum);
        }
        return res;
    }
    rgsw& get_rgsw(size_t row, size_t col) const {
        return get(row, col);
    }
    size_t n_rlwes() const {
        return get_rgsw(0, 0).size();
    }
    size_t n_polys() const {
        return get_rgsw(0, 0).n_polys();
    }
    size_t n_coeffs() const {
        return get_rgsw(0, 0).n_coeffs();
    }

};


// TODO use this instead of vector_i128
class veri_vec_scalar : public array1d<i128, veri_vec_scalar> {
private:
public:
    veri_vec_scalar() : array1d<i128, veri_vec_scalar>() {}
    veri_vec_scalar(size_t N) : array1d<i128, veri_vec_scalar>(N) {}
    ~veri_vec_scalar() {}

    rgsw_vec operator*(const rgsw_mat& other) const {
        size_t rows = other.n_rows();
        size_t cols = other.n_cols();
        assert(rows == size());
        size_t n_rlwes = other.n_rlwes();
        size_t n_polys = other.n_polys();
        size_t n_coeffs = other.n_coeffs();
        rgsw_vec res(cols, n_rlwes, n_polys, n_coeffs);
        for (size_t j = 0; j < cols; j++) {
            rgsw sum(n_rlwes, n_polys, n_coeffs);
            for (size_t i = 0; i < rows; i++) {
                // scalar mult member function from rgsw_vec
                auto val = other.get(i, j) * get(i);
                sum = sum + val;
            }
            res.set(j, sum);
        }
        return res;
    }
    rlwe operator*(const rlwe_vec& other) const {
        size_t n = other.size();
        assert(size() == n);
        size_t n_polys = other.get(0).size();
        size_t n_coeffs = other.get(0).get(0).size();
        rlwe sum(n_polys, n_coeffs);

        for (size_t i = 0; i < n; i++) {
            auto val = other.get(i) * get(i);
            sum = sum + val;
        }
        return sum;
    }
    using array1d<i128, veri_vec_scalar>::pow;
    rlwe pow(const rlwe_vec& other) const {
        size_t n = other.size();
        assert(size() == n);
        size_t n_polys = other.get(0).size();
        size_t n_coeffs = other.get(0).get(0).size();
        rlwe product(n_polys, n_coeffs);
        product.set_coeffs_to_one();
        for (size_t i = 0; i < n; i++) {
            auto val = other.get(i).pow(get(i));
            product = product.group_mult(val);
        }
        return product;
    }
};

i128 random_i128(i128 from, i128 to_inclusive) {
    // random int in range [from, to_inclusive]
    // TODO docs suggest the above inclusive range. Double check
    std::random_device rd;  // Non-deterministic seed
    std::mt19937 gen(rd()); // Mersenne Twister engine
    std::uniform_int_distribution<i128> distrib(from, to_inclusive); // Range: 0 to 100
    i128 random_number = distrib(gen);
    return random_number;
}


struct Params {
private:
    // ======== Controller matrices ========
    std::vector<std::vector<double>> F_ = {
        {2, 0, 0, 0, 0},
        {0, -1, 0, 0, 0},
        {0, 0, 1, 0, 0},
        {0, 0, 0, 0, 0},
        {0, 0, 0, 0, 0}
    };

    std::vector<std::vector<double>> G = {
        {0.0816, 0.0047, 1.6504, -0.0931, 0.4047},
        {-1.4165, -0.3163, -0.4329, 0.1405, 0.8263},
        {-1.4979, -0.2089, -0.6394, 0.3682, 0.7396},
        {0.0459, 0.0152, 1.1004, -0.1187, 0.6563},
        {0.0020, 0.0931, 0.0302, -0.0035, 0.0177}
    };

    std::vector<std::vector<double>> R = {
        {-3.5321, 23.1563},
        {-0.5080, -2.3350},
        {2.5496, 0.9680},
        {0.0436, -1.1227},
        {-0.7560, 0.7144}
    };

    std::vector<std::vector<double>> H = {
        {0.0399, -0.3269, -0.2171, 0.0165, -0.0655},
        {-0.0615, -0.0492, -0.0322, -0.0077, -0.0084}
    };

    // ======== Scale up G, R, and H to integers ========
    array2d<i128> scalar_mat_mult(i128 scalar, std::vector<std::vector<double>>& mat, i128 q, size_t rows, size_t cols) {
        array2d<i128> result(rows, cols);
        for (size_t i = 0; i < rows; ++i) {
            for (size_t j = 0; j < cols; ++j) {
                result.set(i, j, mod_(static_cast<i128>(std::round(scalar * mat[i][j])), q));
            }
        }
        return result;
    }

    // TODO update to array1d<i128> (add as another class?)
    vector_i128 scalar_vec_mult(i128 scalar, vector_double& vec, i128 q) {
        // TODO update to array1d<i128> (add as another class?)
        vector_i128 result(vec.size());
        assert(result.capacity() == result.size());
        for (size_t i = 0; i < vec.size(); ++i) {
            result.at(i) = mod_(static_cast<i128>(std::round(scalar * vec[i])), q);
        }
        return result;
    }

    void generate_field_and_group_params() {
        // TODO use NTL and add funcs from Python implementation
        p = GROUP_MODULUS;
        q = FIELD_MODULUS;
        g = GENERATOR;
    }

public:
    Params() :
        N(N_),
        p(GROUP_MODULUS),
        q(FIELD_MODULUS),
        g(GENERATOR),
        s(10001),
        L(10001),
        r(10001),
        iter_(10),
        F(F_.size(), F_.at(0).size()),
        G_bar(scalar_mat_mult(s, G, q, G.size(), G.at(0).size())),
        R_bar(scalar_mat_mult(s, R, q, R.size(), R.at(0).size())),
        H_bar(scalar_mat_mult(s, H, q, H.size(), H.at(0).size()))
    {
        // TODO update or remove
        assert(A.size() == B.size());
        assert(A.size() == F_.size());
        // generate_field_and_group_params();
        // TODO move into initialiser list
        x_cont_init_scaled = scalar_vec_mult(r * s * L, x_cont_init, q);
        // Construct F from F_
        for (size_t i = 0; i < F_.size(); i++) {
            for (size_t j = 0; j < F_.at(0).size(); j++) {
                F.set(i, j, F_.at(i).at(j));
            }
        }
    }
    i128 N, p, q, g, s, L, r, iter_, x_dim, y_dim, u_dim;
    array2d<i128> F, G_bar, R_bar, H_bar;
    // TODO update to array1d<i128> (add as another class?)
    vector_i128 x_cont_init_scaled;

    // ======== Plant matrices ========
    std::vector<std::vector<double>> A = {
        {1.0, 0.0020, 0.0663, 0.0047, 0.0076},
        {0.0, 1.0077, 2.0328, -0.5496, -0.0591},
        {0.0, 0.0478, 0.9850, -0.0205, -0.0092},
        {0.0, 0.0, 0.0, 0.3679, 0.0},
        {0.0, 0.0, 0.0, 0.0, 0.3679}
    };
    std::vector<std::vector<double>> B = {
        {0.0029, 0.0045},
        {-0.3178, -0.0323},
        {-0.0086, -0.0051},
        {0.6321, 0},
        {0, 0.6321}
    };
    std::vector<std::vector<double>> C = {
        {0, 1, 0, 0, 0},
        {0, -0.2680, 47.7600, -4.5600, 4.4500},
        {1, 0, 0, 0, 0},
        {0, 0, 0, 1, 0},
        {0, 0, 0, 0, 1}
    };




    // ======== Plant and Controller initial state ========
    std::vector<double> x_plant_init = {
        1,
        -1,
        0,
        0.7,
        2,
    };
    std::vector<double> x_cont_init = {
        -0.001,
        0.013,
        0.2,
        -0.02,
        0,
    };
    void print() {
        std::cout << "*** START Params.print ***\n";
        // print these: p, q, g, s, L, r, iter_;
        std::cout << "N: " << print_to_string_i128(N) << "\n";
        std::cout << "p: " << print_to_string_i128(p) << "\n";
        std::cout << "q: " << print_to_string_i128(q) << "\n";
        std::cout << "g: " << print_to_string_i128(g) << "\n";
        std::cout << "s: " << print_to_string_i128(s) << "\n";
        std::cout << "L: " << print_to_string_i128(L) << "\n";
        std::cout << "r: " << print_to_string_i128(r) << "\n";
        std::cout << "iter: " << print_to_string_i128(iter_) << "\n";
        std::cout << "F:\n";
        F.print();
        std::cout << "G_bar:\n";
        G_bar.print();
        std::cout << "R_bar:\n";
        R_bar.print();
        std::cout << "H_bar:\n";
        H_bar.print();
        std::cout << "x_cont_init_scaled:\n";
        print_vector_i128(x_cont_init_scaled);
        std::cout << "*** END Params.print ***\n\n\n";
    }
    // TODO pass in gen as a param
    vector_i128 sample_knowledge_exponents(i128 from, i128 to_inclusive) {
        i128 N = 6;
        vector_i128 res(N);
        for (size_t i = 0; i < N; i++) {
            // TODO update range
            res.at(i) = random_i128(from, to_inclusive);
            DEBUG1(res.at(i) = 1;)
        }
        // res <- {alpha_0, alpha_1, gamma_0, gamma_1, rho_0, rho_1}
        return res;
    }

    // TODO pass in gen as a param
    std::vector<vector_i128> sample_verification_vectors(i128 m, i128 n, i128 from, i128 to_inclusive) {
        i128 N = 3;
        std::vector<vector_i128> res(N);
        for (size_t i = 0; i < N - 1; i++) {
            res.at(i) = vector_i128(n);
        }
        res.at(N - 1) = vector_i128(m);

        for (auto& x : res) {
            for (size_t i = 0; i < x.size(); i++) {
                x.at(i) = random_i128(from, to_inclusive);
                DEBUG1(x.at(i) = 1;)
            }
        }
        // res <- {r_0, r_1, s}
        return res;
    }
};

struct eval_key {
    std::vector<hashed_a_poly> gr;
    hashed_a_rgsw_vec grFr;
    hashed_a_rgsw_vec gsHr;
    std::vector<hashed_a_poly> gr_rho;
    hashed_a_rgsw_vec grFr_alpha;
    hashed_a_rgsw_vec gsHr_gamma;
};

struct veri_key {
    vector_i128 s;
    hashed_rgsw_vec rG_0;
    hashed_rgsw_vec rG_1;
    hashed_rgsw_vec rR_0;
    hashed_rgsw_vec rR_1;
    i128 rho_0;
    i128 rho_1;
    i128 alpha_0;
    i128 alpha_1;
    i128 gamma_0;
    i128 gamma_1;
};

struct Proof {
    hashed_rlwe grx_;
    hashed_rlwe grFrx;
    hashed_rlwe gsHrx;

    hashed_rlwe gr_rho_x_;
    hashed_rlwe grFr_alpha_x;
    hashed_rlwe gsHr_gamma_x;
    hashed_rlwe g_1;
};