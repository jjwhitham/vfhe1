// TODO replace vectori128 by an array1d i128 subclass which:
    // TODO implements scalar mult (should be covered by base class anyway...)
// TODO use unsigned ints so that mod runs quicker
#pragma once

#include <iostream>
#include <cassert>
#include <random>
#include "omp.h"
#include "shared.h"
#include "ntt.h"
#include "gmpxx.h"

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

#ifdef ASSERT_ON
#  define ASSERT(x) assert(x)
#else
#  define ASSERT(x)
#endif

#ifdef CHECK_ON
#  define CHECK(x)
#else
#  define CHECK(x)
#endif

#ifndef N_THREADS
#  define N_THREADS 1
#endif


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
        for (size_t i = 0; i < size_; i++) arr[i] = other.arr[i];
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
        CHECK(check_index_bounds(n);)
        // if (!disable_value_check)
        CHECK(check_value_bounds(arr[n]);)
        return arr[n];
    }
    void check_value_bounds(const T& val) const {
        if constexpr (std::is_same_v<T, i128>) {
            i128 min_val = 0;
            i128 max_val = GROUP_MODULUS - 1;
            if (val < min_val || val > max_val) {
            // if (val > max_val) {
                throw std::out_of_range(
                    "(array1d) Value out of range: " + i128str(val)
                );
            }
        }
    }
    // void set(int n, T val, bool disable_value_check = false) { // FIXME should be T& ?
    void set(int n, const T& val) {
        CHECK(check_index_bounds(n);)
        // if (!disable_value_check)
        CHECK(check_value_bounds(val);)
        arr[n] = val;
    }
    size_t size() const {
        return size_;
    }
    T& operator[](int i) const {
        return get(i);
    }
    // Copy assignment operator
    array1d& operator=(const array1d& other) {
        if (this != &other) {
            delete[] arr;
            size_ = other.size_;
            arr = new T[size_];
            for (size_t i = 0; i < size_; i++) arr[i] = other.arr[i];
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
            T val = get(i) * scalar;
            if constexpr (std::is_same_v<T, i128>)
                val = mod_(val, FIELD_MODULUS);
            res.set(i, val);
        }
        return res;
    }
    mpz group_mult_(mpz a, mpz b) const {
        return mod_(a * b, GROUP_MODULUS);
    }
    Derived group_mult(const Derived& other) const {
        size_t N = size();
        Derived res(N);
        for (size_t i = 0; i < N; i++) {
            T val;
            if constexpr (std::is_same_v<T, mpz>) {
                val = group_mult_(get(i), other.get(i));
            } else {
                val = get(i).group_mult(other.get(i));
            }
            res.set(i, val);
        }
        return res;
    }
    // TODO might need to change i128 to either mpz or U (if i128 & mpz required)
    Derived group_mult(const i128& scalar) const {
        size_t N = size();
        Derived res(N);
        for (size_t i = 0; i < N; i++) {
            T val;
            if constexpr (std::is_same_v<T, i128>) {
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
            T val = get(i) + other.get(i);
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
            T val;
            if constexpr (std::is_same_v<T, i128>)
                val = mod_sub(get(i), other.get(i));
            else
                val = get(i) - other.get(i);
            // constexpr bool is128 = (std::is_same_v<T, i128>);
            // T val = is128 ? mod_sub(get(i), other.get(i)) : get(i) - other.get(i);
            neg_other.set(i, val);
        }
        return neg_other;
    }
    // raise self to other
    Derived pow(const Derived& other) const {
        size_t N = size();
        Derived res(N);
        if constexpr (std::is_same_v<T, mpz>) {
            for (size_t i = 0; i < N; i++) {
                auto val = pow_(get(i), other.get(i), GROUP_MODULUS);
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
    Derived pow(const mpz scalar) const {
        size_t N = size();
        Derived res(N);
        if constexpr (std::is_same_v<T, mpz>) {
            for (size_t i = 0; i < N; i++) {
                auto val = pow_(scalar, get(i), GROUP_MODULUS);
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
        if constexpr (std::is_same_v<T, mpz>) {
            for (size_t i = 0; i < N; i++) {
                auto val = pow_(GENERATOR, get(i), GROUP_MODULUS);
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
    // // base case binary modular exponentiation
    // i128 pow_(i128 base, i128 power) const {
    //     // power = mod_(power, FIELD_MODULUS);
    //     // base = mod_(base, GROUP_MODULUS);
    //     i128 result = 1;
    //     while (power > 0) {
    //         bool is_power_odd = (power % 2) == 1;
    //         if (is_power_odd)
    //             result = (result * base) % GROUP_MODULUS;
    //         power >>= 1;
    //         base = (base * base) % GROUP_MODULUS;
    //     }
    //     return result;
    // }

    // // base case binary modular exponentiation
    // mpz_class pow_(mpz_class base, mpz_class power) const {
    //     mpz result = 1;
    //     // TODO should we just mutate base?
    //     mpz_powm(result.get_mpz_t(), base.get_mpz_t(), power.get_mpz_t(), GROUP_MODULUS.get_mpz_t());
    //     return result;
    // }


    T* begin() { return arr; }
    T* end() { return arr + size_; }
    const T* begin() const { return arr; }
    const T* end() const { return arr + size_; }

    void print() const {
        if constexpr (std::is_same_v<T, i128>) {
            std::cout << "{";
            for (size_t i = 0; i < size() - 1; i++) {
                std::cout << i128str(get(i));
                std::cout << ", ";
            }
            std::cout << i128str(get(size() - 1)) << "}";
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
            // i128 min_val = -1;
            // min_val <<= 63;
            // i128 max_val = (1UL << 63) - 1;
            // if (val < min_val || val > max_val) {
            if (val >= FIELD_MODULUS) {
                throw std::out_of_range(
                    "(array2d) Value out of range: " + i128str(val)
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
    void print_i128_old() const {
        for (size_t i = 0; i < rows_; i++) {
            std::cout << "{";
            for (size_t j = 0; j < cols_ - 1; j++) {
                std::cout << i128str(arr[i][j]) << ", ";
            }
            std::cout << i128str(arr[i][cols_ - 1]) << "}\n";
        }
    }
    void print_i128() const {
        for (size_t i = 0; i < rows_; i++) {
            std::cout << "{";
            for (size_t j = 0; j < cols_ - 1; j++) {
                std::cout << print_to_string_mpz(arr[i][j]) << ", ";
            }
            std::cout << print_to_string_mpz(arr[i][cols_ - 1]) << "}\n";
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
        if (N == 2 * N_) // FIXME do better
            isNTT = true;
    }
    void set_isNTT(bool isNTT_) {
        assert(isNTT != isNTT_);
        isNTT = isNTT_;
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


    poly convolve(const poly& other) const {
        // FIXME n_coeffs differs when first poly is in NTT form
        ASSERT(n_coeffs() == 2 * other.n_coeffs());
        // NOTE rgsw mats should already be converted to NTT form
        assert(isNTT);

        size_t n = POLY_SIZE;
        poly a(2 * n);
        for (size_t i = 0; i < 2 * n; i++)
            a.set(i, get(i));

        poly b(2 * n);
        for (size_t i = 0; i < n; i++)
            b.set(i, other.get(i));
        arr_u128 psi_pows = get_rou_pows(TWO_ROU);
        ntt_iter(b, psi_pows);

        for (size_t i = 0; i < 2 * POLY_SIZE; i++) {
            a.set(i, (a.get(i) * b.get(i)) % FIELD_MOD);
        }
        a.isNTT = true;

        return a;
    }
    poly nega_ntt(const poly& other) const {
        assert(n_coeffs() == other.n_coeffs());
        size_t n = POLY_SIZE;
        // NOTE constructor zero-initialises
        poly a(n);
        poly b(n);
        for (size_t i = 0; i < n; i++) {
            a.set(i, get(i));
            b.set(i, other.get(i));
        }

        u128 INV_2ROU = pow_constexpr(TWO_ROU, FIELD_MOD - 2, FIELD_MOD);
        u128 INV_N = pow_constexpr(POLY_SIZE, FIELD_MOD - 2, FIELD_MOD);
        arr_u128 psi_pows = get_rou_pows(TWO_ROU);
        arr_u128 psi_inv_pows = get_rou_pows(INV_2ROU);
        ntt_iter1(a, psi_pows);
        ntt_iter1(b, psi_pows);
        for (size_t i = 0; i < POLY_SIZE; i++) {
            a.set(i, (a.get(i) * b.get(i)) % FIELD_MOD);
        }
        intt_iter1(a, psi_inv_pows, INV_N);
        return a;
    }

    poly conv_to_ntt() {
        ASSERT(isNTT == false);
        size_t n = n_coeffs();
        // NOTE constructor zero-initialises
        poly a(2 * n);
        for (size_t i = 0; i < n; i++) {
            a.set(i, get(i));
        }
        arr_u128 psi_pows = get_rou_pows(TWO_ROU);
        ntt_iter(a, psi_pows);
        a.isNTT = true;
        return a;
    }
    poly conv_to_coeff() {
        ASSERT(isNTT == true);
        size_t n = n_coeffs();
        // NOTE constructor zero-initialises
        poly a(n);
        for (size_t i = 0; i < n; i++) {
            a.set(i, get(i));
        }
        u128 INV_2ROU = pow_constexpr(TWO_ROU, FIELD_MOD - 2, FIELD_MOD);
        arr_u128 psi_inv_pows = get_rou_pows(INV_2ROU);
        u128 INV_2N = pow_constexpr(2 * POLY_SIZE, FIELD_MOD - 2, FIELD_MOD);
        intt_iter(a, psi_inv_pows, INV_2N);
        a.isNTT = false;
        assert(a.get(n - 1) == 0);
        return a;
    }

    using array1d<i128, poly>::operator*;
    poly operator*(const poly& other) const {
        ASSERT(n_coeffs() == other.n_coeffs());
        return nega_ntt(other);
    }
    // FIXME - remove dead code
    poly conv_to_nega_(poly& conv) const {
        // std::cout << conv.size() << "\n";
        size_t n = conv.size() / 2;
        assert(conv.get(2 * n - 1) == 0);
        poly res(n);
        for (size_t i = 0; i < n - 1; i++) {
            // i128 val = mod_(conv.get(i) - conv.get(i + n), FIELD_MODULUS);
            i128 val = mod_sub(conv.get(i), conv.get(i + n));
            res.set(i, val);
        }
        res.set(n - 1, conv.get(n - 1));
        return res;
    }
    auto conv_to_nega(size_t N) const {
        TIMING(timing.calls_conv_to_nega += 1;)
        // ASSERT(n_coeffs() == 2 * N - 1);
        size_t conv_degree = n_coeffs() - 1;
        // HACK can't return *this after defining array1d move semantics
        // ASSERT(conv_degree >= N);
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

class hashed_a_poly : public array1d<mpz, hashed_a_poly> {
private:
public:
    hashed_a_poly() : array1d<mpz, hashed_a_poly>() {}
    hashed_a_poly(size_t n_hashed_a_coeffs) : array1d<mpz, hashed_a_poly>(n_hashed_a_coeffs) {
        for (size_t i = 0; i < n_hashed_a_coeffs; i++)
            set(i, 0);
    }
    auto& get_hashed_a_coeff(size_t n) const {
        return get(n);
    }
    size_t n_hashed_a_coeffs() const {
        return size();
    }
    mpz get_hash_sec(const poly& other) const {
        TIMING(auto start = std::chrono::high_resolution_clock::now();)
        TIMING(timing.calls_get_hash_sec += 1;)
        // HACK get around gr having hashed_a_polys of length 2n-1, but x_nega_ only length n
        // ASSERT(size() == other.size() || size() == (2 * other.size() - 1));
        mpz result = 1;
        std::array<mpz, N_THREADS> partials;
        for (size_t t = 0; t < N_THREADS; t++) partials[t] = 1;
        #pragma omp parallel for num_threads(N_THREADS)
        for (size_t i = 0; i < other.size(); i++) {
            int tid = omp_get_thread_num();
            partials[tid] = group_mult_(partials[tid], pow_(get(i), other.get(i), GROUP_MODULUS));
            // mpz_powm(res.get_mpz_t(), res.get_mpz_t(), exp.get_mpz_t(), GROUP_MODULUS.get_mpz_t());
            // partials[tid] = group_mult_(partials[tid], mpz_powm(get(i)., other.get(i));
        }
        for (size_t t = 0; t < N_THREADS; t++) {
            result = group_mult_(result, partials[t]);
        }
        TIMING(auto end = std::chrono::high_resolution_clock::now();)
        TIMING(timing.get_hash_sec += end - start;)
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
        ASSERT(n_hashed_a_polys == N_POLYS_IN_RLWE);
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
};

class hashed_rlwe : public array1d<i128, hashed_rlwe> {
private:
public:
    // NOTE hashed_rlwe always has two hashed_polys
    hashed_rlwe() : array1d<i128, hashed_rlwe>() {}
    hashed_rlwe(size_t n_hashed_polys) : array1d<i128, hashed_rlwe>(n_hashed_polys) {
        ASSERT(n_hashed_polys == N_POLYS_IN_RLWE);
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
        ASSERT(n_polys == 2);
    }
    rlwe(size_t n_polys, size_t n_coeffs) : array1d<poly, rlwe>(n_polys) {
        ASSERT(n_polys == 2);
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
    rlwe_decomp decompose(const u128& v, const u32& d, const u32& power) const {
        rlwe_decomp polys(2 * d, n_coeffs());
        for (size_t k = 0; k < N_POLYS_IN_RLWE; k++) {
            poly pol = get_poly(k);
            for (size_t i = 0; i < n_coeffs(); i++) {
                for (size_t j = 0; j < (size_t)d; j++) {
                    i128 decomped_coeff = (v - 1) & (pol.get_coeff(i) >> (power * j));
                    ASSERT(decomped_coeff < v);
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
    void conv_to_ntt() {
        for (size_t i = 0; i < N_POLYS_IN_RLWE; i++) {
            set(i, get_poly(i).conv_to_ntt());
        }
    }
    void conv_to_coeff() {
        for (size_t i = 0; i < N_POLYS_IN_RLWE; i++) {
            set(i, get_poly(i).conv_to_coeff());
        }
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
            ASSERT(n_hashed_polys == 2);
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
    auto decompose(i128 v, u32 d, u32 power) const {
        rlwe_decomp_vec decomps(size(), 2 * d, n_coeffs());
        for (size_t i = 0; i < size(); i++)
            decomps.set(i, get_rlwe(i).decompose(v, d, power));
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
    void conv_to_coeff() const {
        for (size_t i = 0; i < n_rlwes(); i++) {
            get_rlwe(i).conv_to_coeff();
        }
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
        // size_t n_polys_ = n_hashed_a_polys();
        ASSERT(n_hashed_a_polys() == N_POLYS_IN_RLWE);
        ASSERT(n_hashed_a_rlwes() == other.size());
        hashed_rlwe_vec res_vec(n_hashed_a_rlwes(), n_hashed_a_polys());
        // #pragma omp parallel for num_threads(N_THREADS)
        for (size_t i = 0; i < n_hashed_a_rlwes(); i++) {
            hashed_rlwe& res = res_vec.get(i);
            for (size_t j = 0; j < n_hashed_a_polys(); j++) {
                auto val = get(i).get(j).get_hash_sec(other.get(i));
                res.set(j, val);
            }
        }

        hashed_rlwe res(n_hashed_a_polys());
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
        ASSERT(n_hashed_polys() == N_POLYS_IN_RLWE);
        ASSERT(n_hashed_rlwes() == other.n_hashed_polys());
        hashed_rlwe_vec res_vec(n_hashed_rlwes(), n_hashed_polys());
        // #pragma omp parallel for num_threads(2)
        for (size_t i = 0; i < n_hashed_rlwes(); i++) {
            // FIXME a bit hacky - obj.pow(scalar) is scalar^obj, do we need obj^scalar???
            hashed_rlwe& res = res_vec.get(i);
            for (size_t j = 0; j < n_hashed_polys(); j++) {
                auto val = pow_(get(i).get(j), (other.get(i)), GROUP_MODULUS);
                res.set(j, val);
            }
        }

        hashed_rlwe res(n_hashed_polys());
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
        ASSERT(n_rlwes() == other.n_polys());
        rlwe res(n_polys(), n_coeffs());
        for (size_t i = 0; i < n_rlwes(); i++) {
            for (size_t j = 0; j < n_polys(); j++) {
                poly& p = res.get(j);
                auto val = get(i).get(j) * other.get(i);
                p = p + val;
                // res.set(j, res.get(j) + val);
            }
        }
        return res;
    }

    rlwe convolve(const rlwe_decomp& other) const {
        ASSERT(n_rlwes() == other.n_polys());
        rlwe res(n_polys(), 2 * POLY_SIZE); // FIXME
        // TODO need thread array to accumulate sum
        // #pragma omp parallel for num_threads(N_THREADS)
        for (size_t i = 0; i < n_rlwes(); i++) {
            for (size_t j = 0; j < n_polys(); j++) {
                poly& p = res.get(j);
                auto val = get(i).get(j).convolve(other.get(i));
                p = p + val;
            }
        }
        return res;
    }

    void conv_to_ntt() {
        for (size_t i = 0 ; i < n_rlwes(); i++)
            get_rlwe(i).conv_to_ntt();
    }

    using array1d<rlwe, rgsw>::pow;
    rlwe pow(const rlwe_decomp& other) const {
        ASSERT(n_rlwes() == other.n_polys());
        rlwe_vec res_vec(n_rlwes(), n_polys(), n_coeffs());
        // #pragma omp parallel for num_threads(N_THREADS)
        for (size_t i = 0; i < n_rlwes(); i++) {
            rlwe& res = res_vec.get(i);
            for (size_t j = 0; j < n_polys(); j++) {
                auto val = get(i).get(j).pow(other.get(i));
                res.set(j, val);
            }
        }

        rlwe res(n_polys(), n_coeffs());
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
        ASSERT(n_hashed_polys == 2);
        return n_hashed_polys;
    }
    using array1d<hashed_rgsw, hashed_rgsw_vec>::pow;
    hashed_rlwe pow(const hashed_rlwe_decomp_vec& other) const {
        ASSERT(n_hashed_rgsws() == other.n_hashed_rlwe_decomps());
        ASSERT(n_hashed_polys() == 2);
        hashed_rlwe_vec res_vec(n_hashed_rgsws(), n_hashed_polys());
        hashed_rlwe res(n_hashed_polys());
        res.set_coeffs_to_one();
        // #pragma omp parallel for num_threads(N_THREADS)
        for (size_t i = 0; i < n_hashed_rgsws(); i++) {
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
        ASSERT(n_hashed_a_polys == 2);
        return n_hashed_a_polys;
    }
    // using array1d<hashed_a_rgsw, hashed_a_rgsw_vec>::pow;
    hashed_rlwe get_hash_sec(const rlwe_decomp_vec& other) const {
        ASSERT(n_hashed_a_rgsws() == other.n_rlwe_decomps());
        ASSERT(n_hashed_a_polys() == 2);
        hashed_rlwe_vec res_vec(n_hashed_a_rgsws(), n_hashed_a_polys());
        hashed_rlwe res(n_hashed_a_polys());
        res.set_coeffs_to_one();
        // #pragma omp parallel for num_threads(N_THREADS)
        for (size_t i = 0; i < n_hashed_a_rgsws(); i++) {
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
        ASSERT(n_rgsws() == other.n_rlwe_decomps());
        rlwe sum(n_polys(), n_coeffs());
        for (size_t i = 0; i < n_rgsws(); i++) {
            auto val = get(i) * other.get(i);
            sum = sum + val;
        }
        return sum;
    }
    using array1d<rgsw, rgsw_vec>::pow;
    rlwe pow(const rlwe_decomp_vec& other) const {
        ASSERT(n_rgsws() == other.n_rlwe_decomps());
        rlwe_vec res_vec(n_rgsws(), n_polys(), n_coeffs());
        rlwe res(n_polys(), n_coeffs());
        res.set_coeffs_to_one();
        // #pragma omp parallel for num_threads(N_THREADS)
        for (size_t i = 0; i < n_rgsws(); i++) {
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
        ASSERT(n_cols() == other.size());
        ASSERT(n_coeffs() == other.n_coeffs());

        rlwe_vec res(n_rows(), n_polys(), n_coeffs());
        for (size_t i = 0; i < n_rows(); i++) {
            rlwe sum(n_polys(), n_coeffs());
            for (size_t j = 0; j < n_cols(); j++) {
                auto val = get(i, j) * other.get(j);
                sum = sum + val;
            }
            res.set(i, sum);
        }
        return res;
    }
    rlwe_vec convolve(const rlwe_decomp_vec& other) const {
        ASSERT(n_cols() == other.size());
        ASSERT(POLY_SIZE == other.n_coeffs()); // FIXME

        rlwe_vec res(n_rows(), n_polys(), 2 * n_coeffs());
        for (size_t i = 0; i < n_rows(); i++) {
            rlwe sum(n_polys(), 2 * POLY_SIZE); // FIXME
            for (size_t j = 0; j < n_cols(); j++) {
                auto val = get(i, j).convolve(other.get(j));
                sum = sum + val;
            }
            res.set(i, sum);
        }
        return res;
    }
    rlwe_vec pow(const rlwe_decomp_vec& other) const {
        ASSERT(n_cols() == other.size());

        rgsw& rg = get(0, 0);
        rlwe_vec res(n_rows(), rg.n_polys(), rg.n_coeffs());
        for (size_t i = 0; i < n_rows(); i++) {
            rlwe sum(rg.n_polys(), rg.n_coeffs());
            sum.set_coeffs_to_one();
            for (size_t j = 0; j < n_cols(); j++) {
                auto val = get(i, j).pow(other.get(j));
                sum = sum.group_mult(val);
            }
            res.set(i, sum);
        }
        return res;
    }
    void conv_to_ntt() {
        for (size_t row = 0; row < n_rows(); row++) {
            for (size_t col = 0; col < n_cols(); col++) {
                get_rgsw(row, col).conv_to_ntt();
            }
        }
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
        const rgsw_mat& o = other;
        ASSERT(o.n_rows() == size());
        rgsw_vec res(o.n_cols(), o.n_rlwes(), o.n_polys(), o.n_coeffs());
        for (size_t j = 0; j < o.n_cols(); j++) {
            rgsw sum(o.n_rlwes(), o.n_polys(), o.n_coeffs());
            for (size_t i = 0; i < o.n_rows(); i++) {
                auto val = o.get(i, j) * get(i);
                sum = sum + val;
            }
            res.set(j, sum);
        }
        return res;
    }
    rlwe operator*(const rlwe_vec& other) const {
        size_t n = other.size();
        ASSERT(size() == n);
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
        ASSERT(size() == n);
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
    // HACK
    long int from_ = mpz_get_si(from.get_mpz_t());
    long int to_inclusive_ = mpz_get_si(to_inclusive.get_mpz_t());
    std::random_device rd;  // Non-deterministic seed
    std::mt19937 gen(rd()); // Mersenne Twister engine
    std::uniform_int_distribution<unsigned long> distrib(from_, to_inclusive_); // Range: 0 to 100
    mpz random_number = distrib(gen);
    return random_number;
}

u32 get_decomp_power() {
    // v - power of 2, s.t. v^{d-1} < q < v^d
    double log2mpz = log2_mpz(FIELD_MODULUS);
    size_t q_bits = static_cast<size_t>(log2mpz) + 1;
    u32 remainder = q_bits % N_DECOMP;
    u32 power = remainder ? q_bits / N_DECOMP + 1 : q_bits / N_DECOMP;
    assert(FIELD_MODULUS <= mpz(1) << (power * N_DECOMP));
    assert(FIELD_MODULUS > mpz(1) << (power * (N_DECOMP - 1)));
    return power;
}

struct Params {
private:
    // ======== Controller matrices ========
    matrix_double F_ = {
        {2, 0, 0, 0, 0},
        {0, -1, 0, 0, 0},
        {0, 0, 1, 0, 0},
        {0, 0, 0, 0, 0},
        {0, 0, 0, 0, 0}
    };

    matrix_double G = {
        {0.0816, 0.0047, 1.6504, -0.0931, 0.4047},
        {-1.4165, -0.3163, -0.4329, 0.1405, 0.8263},
        {-1.4979, -0.2089, -0.6394, 0.3682, 0.7396},
        {0.0459, 0.0152, 1.1004, -0.1187, 0.6563},
        {0.0020, 0.0931, 0.0302, -0.0035, 0.0177}
    };

    matrix_double R = {
        {-3.5321, 23.1563},
        {-0.5080, -2.3350},
        {2.5496, 0.9680},
        {0.0436, -1.1227},
        {-0.7560, 0.7144}
    };

    matrix_double H = {
        {0.0399, -0.3269, -0.2171, 0.0165, -0.0655},
        {-0.0615, -0.0492, -0.0322, -0.0077, -0.0084}
    };

    // FIXME needs to be mapped from negative to [0,q)
    // ======== Scale up G, R, and H to integers ========
    array2d<i128> scalar_mat_mult(double scalar, matrix_double& mat, i128 q, size_t rows, size_t cols) {
        array2d<i128> result(rows, cols);
        for (size_t i = 0; i < rows; i++) {
            for (size_t j = 0; j < cols; j++) {
                mpf_class rounded = scalar * mat[i][j];
                rounded = mpf_round(rounded);
                i128 modded;
                if (rounded < 0)
                    modded = q + rounded;
                else
                    modded = rounded;
                result.set(i, j, modded);
            }
        }
        return result;
    }

    // FIXME extract lambdas from vfhe.cpp and remove this double up
    // TODO update to array1d<i128> (add as another class?)
    vector_i128 scalar_vec_mult(double scalar, vector_double& vec, i128 q) {
        // TODO update to array1d<i128> (add as another class?)
        vector_i128 result(vec.size());
        ASSERT(result.capacity() == result.size());
        for (size_t i = 0; i < vec.size(); i++) {
            mpf_class rounded = mpf_round(scalar * vec[i]);
            i128 modded;
            if (rounded < 0)
                modded = rounded + q;
            else
                modded = rounded;
            result.at(i) = modded;
        }
        return result;
    }

    // void generate_field_and_group_params() {
    //     // TODO use NTL and add funcs from Python implementation
    //     p = GROUP_MODULUS;
    //     q = FIELD_MODULUS;
    //     g = GENERATOR;
    // }

public:
    Params() :
        N(N_),
        iter_(1),
        s(10000.0),
        L(10000.0),
        r(10000.0),
        p(GROUP_MODULUS),
        q(FIELD_MODULUS),
        g(GENERATOR),
        F(F_.size(), F_.at(0).size()),
        G_bar(scalar_mat_mult(s, G, q, G.size(), G.at(0).size())),
        R_bar(scalar_mat_mult(s, R, q, R.size(), R.at(0).size())),
        H_bar(scalar_mat_mult(s, H, q, H.size(), H.at(0).size()))
    {
        // TODO update or remove
        ASSERT(A.size() == B.size());
        ASSERT(A.size() == F_.size());
        // generate_field_and_group_params();
        // TODO move into initialiser list
        x_cont_init_scaled = scalar_vec_mult(r * s * L, x_cont_init, q);
        // Construct F from F_
        for (size_t i = 0; i < F_.size(); i++) {
            for (size_t j = 0; j < F_.at(0).size(); j++) {
                mpf_class val = F_.at(i).at(j);
                i128 val1;
                if (val < 0)
                    val1 = val + FIELD_MODULUS;
                else
                    val1 = val;
                F.set(i, j, val1);
            }
        }
        assert(d >= 1);
    }
    size_t N, iter_;
    double s, L, r;
    i128 p, q, g;

    array2d<i128> F, G_bar, R_bar, H_bar;
    // TODO update to array1d<i128> (add as another class?)
    vector_i128 x_cont_init_scaled;
    u32 d = N_DECOMP;
    u32 power = get_decomp_power();
    // static constexpr i128 v = static_cast<i128>(1) << power;
    u32 v = 1 << power;

    // ======== Plant matrices ========
    matrix_double A = {
        {1.0, 0.0020, 0.0663, 0.0047, 0.0076},
        {0.0, 1.0077, 2.0328, -0.5496, -0.0591},
        {0.0, 0.0478, 0.9850, -0.0205, -0.0092},
        {0.0, 0.0, 0.0, 0.3679, 0.0},
        {0.0, 0.0, 0.0, 0.0, 0.3679}
    };
    matrix_double B = {
        {0.0029, 0.0045},
        {-0.3178, -0.0323},
        {-0.0086, -0.0051},
        {0.6321, 0},
        {0, 0.6321}
    };
    matrix_double C = {
        {0, 1, 0, 0, 0},
        {0, -0.2680, 47.7600, -4.5600, 4.4500},
        {1, 0, 0, 0, 0},
        {0, 0, 0, 1, 0},
        {0, 0, 0, 0, 1}
    };




    // ======== Plant and Controller initial state ========
    vector_double x_plant_init = {
        1,
        -1,
        0,
        0.7,
        2,
    };
    vector_double x_cont_init = {
        -0.001,
        0.013,
        0.2,
        -0.02,
        0,
    };
    void print() {
        std::cout << "*** START Params.print ***\n";
        // print these: p, q, g, s, L, r, iter_;
        std::cout << "N: " << i128str(N) << "\n";
        // std::cout << "p: " << i128str(p) << "\n";
        std::cout << "p: " << print_to_string_mpz(p) << "\n";
        // std::cout << "q: " << i128str(q) << "\n";
        std::cout << "q: " << print_to_string_mpz(q) << "\n";
        // std::cout << "g: " << i128str(g) << "\n";
        std::cout << "g: " << print_to_string_mpz(g) << "\n";
        std::cout << "s: " << i128str(s) << "\n";
        std::cout << "L: " << i128str(L) << "\n";
        std::cout << "r: " << i128str(r) << "\n";
        std::cout << "iter: " << i128str(iter_) << "\n";
        std::cout << "F:\n";
        F.print();
        std::cout << "G_bar:\n";
        G_bar.print();
        std::cout << "R_bar:\n";
        R_bar.print();
        std::cout << "H_bar:\n";
        H_bar.print();
        std::cout << "x_cont_init_scaled:\n";
        // print_vector_i128(x_cont_init_scaled);
        print_vector_mpz(x_cont_init_scaled);
        std::cout << "*** END Params.print ***\n\n\n";
    }
    // TODO pass in gen as a param
    vector_i128 sample_knowledge_exponents(i128 from, i128 to_inclusive) {
        size_t N = 6;
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
    std::vector<vector_i128> sample_verification_vectors(__uint128_t m, __uint128_t n, i128 from, i128 to_inclusive) {
        size_t N = 3;
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
    const matrix_double& get_F_() const {
        return F_;
    }
    const matrix_double& get_G() const {
        return G;
    }
    const matrix_double& get_H() const {
        return H;
    }
    const matrix_double& get_R() const {
        return R;
    }
    const vector_double& get_x_cont_init() const {
        return x_cont_init;
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