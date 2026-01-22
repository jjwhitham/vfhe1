/* XXX
void millerLoopVec(GT& z, const G1 *x, const G2 *y, size_t n);
TODO
    This function is for multi-pairing
        computes prod_{i=0}^{n-1} MillerLoop(x[i], y[i])
        prod_{i=0}^{n-1} e(x[i], y[i]) = finalExp(prod_{i=0}^{n-1} MillerLoop(x[i], y[i]))
*/
#pragma once

#include <iostream>
#include <cassert>
#include <random>
#include <omp.h>
#include "shared.h"
#include "ntt.h"
#include "gmpxx.h"

// #include <boost/stacktrace.hpp>
// #include <boost/exception/all.hpp>
// typedef boost::error_info<struct tag_stacktrace, boost::stacktrace::stacktrace> traced;

#include <typeinfo>
#ifdef __GNUG__
#include <cxxabi.h>
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

#ifdef ASSERT_ON
#  define ASSERT(x) assert(x)
#else
#  define ASSERT(x)
#endif

#ifdef CHECK_ON
#  define CHECK(x) x
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
    array1d(size_t size) : size_{size}, arr{new T[size]} {}
    // Copy constructor
    array1d(const array1d& other) : size_(other.size_), arr(new T[other.size_]) {
        for (size_t i = 0; i < size_; i++)
            arr[i] = other.arr[i];
    }
    // Move constructor
    array1d(array1d&& other) noexcept : size_(other.size_), arr(other.arr) {
        other.size_ = 0;
        other.arr = nullptr;
    }
    ~array1d() {
        delete[] arr;
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
    void check_index_bounds(size_t n) const {
        if (n >= size_) {
            // std::clog << boost::current_exception_diagnostic_information() << std::endl;
            try {
                std::out_of_range e = std::out_of_range(
                    "Index error: accessing arr[" + std::to_string(n) + "]"
                    + " in a " + std::to_string(size_) + " element array."
                );
                // throw boost::enable_error_info(e) << traced(boost::stacktrace::stacktrace());
            } catch (const std::exception& e1) {
                // const boost::stacktrace::stacktrace* st = boost::get_error_info<traced>(e1);
                // if (st) {
                //     std::cerr << *st << '\n';
                // }
            }
        }
    }
    // TODO check ref return is most appropriate
    T& get(size_t n, bool disable_value_check=false) const {
        CHECK(check_index_bounds(n);)
        if (!disable_value_check) {
            CHECK(check_value_bounds(arr[n]);)
        }
        return arr[n];
    }
    void check_value_bounds(const T& val) const {
        if constexpr (std::is_same_v<T, bigz>) {
            static const bigz min_val = 0;
            static const bigz max_val = FIELD_MODULUS - 1;
            if (val < min_val || val > max_val) {
                std::string message = "(array1d) Value out of range: " + std::string(print_to_string_mpz(val)) + "\n";
            #ifdef ASSERT_ON
                throw std::out_of_range(message);
            #else
                std::cout << message;
            #endif
            }
        }
    }
    // void set(int n, T val, bool disable_value_check = false) { // FIXME should be T& ?
    // TODO check that const ref is most appropriate
    void set(int n, const T& val) {
        CHECK(check_index_bounds(n);)
        // if (!disable_value_check)
        CHECK(check_value_bounds(val);)
        arr[n] = val;
    }
    void set(int n, const T&& val) {
        CHECK(check_index_bounds(n);)
        // if (!disable_value_check)
        CHECK(check_value_bounds(val);)
        arr[n] = val;
    }
    size_t size() const {
        return size_;
    }
    T& operator[](int i) const {
        bool disable_value_check = false;
    #ifdef ASSERT_ON
        disable_value_check = true;
    #endif
        return get(i, disable_value_check);
    }
    Derived operator*(const bigz& scalar) const {
        // size_t N = size();
        // Derived res(N);
        // for (size_t i = 0; i < N; i++) {
        //     T val = get(i) * scalar;
        //     if constexpr (std::is_same_v<T, bigz>)
        //         val = mod_(val, FIELD_MODULUS);
        //     res.set(i, val);
        // }
        // return res;

        size_t N = size();
        Derived res(N);
        if constexpr (std::is_same_v<T, bigz>) {
            #pragma omp parallel for schedule(static) num_threads(N_THREADS)
            for (size_t i = 0; i < N; i++) {
                res[i] = mod_((*this)[i] * scalar, FIELD_MODULUS);
            }
        } else {
            for (size_t i = 0; i < N; i++) {
                res[i] = (*this)[i] * scalar;
            }
        }
        return res;
    }
    G1 group_mult_(G1& a, const G1& b) const {
        return a + b;
    }
    /**
     *  @brief Group multiplication.
     *  Multiplies two group elements pointwise and recursively
     *
     *  @param other The other array1d of the same concrete type
     *
     *  @return An array1d of the same type as `this` and `other`
     *
     *  @note Deprecated...
     */
    Derived group_mult(const Derived& other) const {
        size_t N = size();
        Derived res(N);
        for (size_t i = 0; i < N; i++) {
            T val;
            if constexpr (std::is_same_v<T, bigz>) {
                throw std::runtime_error("Derived group_mult(Derived): Received a bigz arg!\n");
            } else if constexpr (std::is_same_v<T, G1>) {
                val = group_mult_(get(i), other.get(i));
            } else {
                val = get(i).group_mult(other.get(i));
            }
            res.set(i, val);
        }
        return res;
    }
    // return a readable name for Derived (demangled on GNU)
    static std::string name() {
        const char* tname = typeid(Derived).name();
    #ifdef __GNUG__
        int status = 0;
        char* dem = abi::__cxa_demangle(tname, nullptr, nullptr, &status);
        std::string res;
        if (status == 0 && dem) {
            res = dem;
            std::free(dem);
        } else {
            res = tname;
        }
        return res;
    #else
        return std::string(tname);
    #endif
    }

    void check() const {
        size_t N = size();
        if constexpr (std::is_same_v<T, bigz>) {
            // std::cout << "below `if constexpr'\n";
            for (size_t i = 0; i < N; i++) {
                // if (i == 0) {
                //     std::cout << name() << ": " << "Checking value bounds\n";
                // }
                T val = arr[i];
                check_value_bounds(val);
            }
        } else {
            for (size_t i = 0; i < N; i++) {
                // if (i == 0) {
                //     std::cout << "hello from " << name() << "\n";
                // }
                arr[i].check();
            }
        }
    }
    void mod() {
        size_t N = size();
        if constexpr (std::is_same_v<T, bigz>) {
            static const bigz min_val = 0;
            static const bigz max_val = FIELD_MODULUS - 1;
            for (size_t i = 0; i < N; i++) {
                bigz& val = arr[i];
                if (val < min_val || val > max_val)
                    val %= FIELD_MODULUS;
                if (val < min_val)
                    val += FIELD_MODULUS;
            }
        } else {
            for (size_t i = 0; i < N; i++) {
                arr[i].mod();
            }
        }
    }
    Derived operator+(const Derived& other) const {
        size_t N = size();
        Derived res(N);
        if constexpr (std::is_same_v<T, bigz>) {
        #pragma omp parallel for schedule(static) num_threads(N_THREADS)
            for (size_t i = 0; i < N; i++) {
                res[i] = mod_((*this)[i] + other[i], FIELD_MODULUS);
            }
        } else {
            for (size_t i = 0; i < N; i++) {
                res[i] = (*this)[i] + other[i];
            }
        }
        return res;
    }
    Derived operator+(Derived&& other) const {
        return operator+(other);
    }
    Derived& operator+=(const Derived& other) {
        size_t N = size();
        if constexpr (std::is_same_v<T, bigz>) {
            #pragma omp parallel for schedule(static) num_threads(N_THREADS)
            for (size_t i = 0; i < N; i++) {
                get(i) += other.get(i);
                    mod_(get(i), FIELD_MODULUS);
            }
        } else {
            for (size_t i = 0; i < N; i++)
                get(i) += other.get(i);
        }
        return static_cast<Derived&>(*this);
    }
    Derived& operator+=(const Derived&& other) {
        size_t N = size();
        if constexpr (std::is_same_v<T, bigz>) {
            #pragma omp parallel for schedule(static) num_threads(N_THREADS)
            for (size_t i = 0; i < N; i++) {
                get(i) += other.get(i);
                    mod_(get(i), FIELD_MODULUS);
            }
        } else {
            for (size_t i = 0; i < N; i++)
                get(i) += other.get(i);
        }
        return static_cast<Derived&>(*this);
    }
    Derived operator-(const Derived& other) const {
        // size_t N = size();
        // Derived neg_other(N);
        // for (size_t i = 0; i < N; i++) {
        //     T val;
        //     if constexpr (std::is_same_v<T, bigz>)
        //         val = mod_sub(get(i), other.get(i));
        //     else
        //         val = get(i) - other.get(i);
        //     neg_other.set(i, val);
        // }
        // return neg_other;
        size_t N = size();
        Derived res(N);
        if constexpr (std::is_same_v<T, bigz>) {
            #pragma omp parallel for schedule(static) num_threads(N_THREADS)
            for (size_t i = 0; i < N; i++) {
                res[i] = mod_sub((*this)[i], other[i]);
            }
        } else {
            for (size_t i = 0; i < N; i++) {
                res[i] = (*this)[i] - other[i];
            }
        }
        return res;
    }

    T* begin() { return arr; }
    T* end() { return arr + size_; }
    const T* begin() const { return arr; }
    const T* end() const { return arr + size_; }

    void print() const {
        if constexpr (std::is_same_v<T, bigz>) {
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
    ~array2d() {
        // delete[] data;
        // delete[] arr;
    } // FIXME
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
        if constexpr (std::is_same_v<T, bigz>) {
            // bigz min_val = -1;
            // min_val <<= 63;
            // bigz max_val = (1UL << 63) - 1;
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
        if constexpr (std::is_same_v<T, bigz>) {
            print_i128();
        } else {
            print_array1d();
        }
    }
};

class poly : public array1d<bigz, poly> {
private:
    bool isNTT = false;
public:
    poly() : array1d<bigz, poly>() {}
    poly(size_t N) : array1d<bigz, poly>(N) {
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
    auto get_hash(const vector_bigz& eval_pows) const {
        // std::cout << "poly::get_hash: size(): " << size() << "\n";
        bigz hash = 0;
        for (size_t i = 0; i < size(); i++) {
            hash += get(i) * eval_pows.at(i);
        }
        mod_(hash, FIELD_MODULUS);
        return hash;
    }
    // XXX See inline definition below class hashed_t_poly (avoids circular definitions)
    auto get_hash_t(const vector_bigz& eval_pows) const;

    poly convolve(const poly& other) const {
        // NOTE rgsw mats should already be converted to NTT form
        assert(isNTT);
        assert(other.isNTT);
        poly a(2 * N_);
        #pragma omp parallel for schedule(static) num_threads(N_THREADS)
        for (size_t i = 0; i < 2 * N_; i++) {
            // XXX = is faster than += here
            a[i] = mod_((*this)[i] * other.get(i), FIELD_MODULUS);
        }
        a.isNTT = true;
        return a;
    }
    // NOTE only used for Enc/Dec
    poly nega_ntt(const poly& other) const {
        assert(other.isNTT == true);
        assert(n_coeffs() == other.n_coeffs());
        poly a{*this}; // NOTE copy constructor

        // TODO should be in ntt.h
        static const arr_u128 psi_pows = get_rou_pows(TWO_ROU);
        ntt_iter1(a, psi_pows);
        #pragma omp parallel for schedule(static) num_threads(N_THREADS)
        for (size_t i = 0; i < N_; i++) {
            // TODO consider using * as pointwise mult and call nega_ntt explicitly
            a[i] = mod_(a[i] * other[i], FIELD_MODULUS);
        }
        static const bigz INV_2ROU = pow_constexpr(TWO_ROU, FIELD_MODULUS - 2, FIELD_MODULUS);
        static const arr_u128 psi_inv_pows = get_rou_pows(INV_2ROU);
        static const bigz INV_N = pow_constexpr(N_, FIELD_MODULUS - 2, FIELD_MODULUS);
        intt_iter1(a, psi_inv_pows, INV_N);
        return a;
    }

    // TODO should mutate
    poly to_eval_form(bool is_conv=true) const {
        ASSERT(isNTT == false);
        size_t n = n_coeffs();
        if (is_conv)
            n *= 2;
        poly a{n};
        // NOTE don't use copy constructor, as conv requires zero-padding to 2n
        for (size_t i = 0; i < n_coeffs(); i++) {
            a.set(i, get(i));
        }
        static const arr_u128 psi_pows = get_rou_pows(TWO_ROU);
        is_conv ? ntt_iter(a, psi_pows) : ntt_iter1(a, psi_pows);
        a.isNTT = true;
        return a;
    }
    poly to_coeff_form(bool is_conv=true) const {
        ASSERT(isNTT == true);
        size_t n = n_coeffs();
        poly a{*this};
        static const bigz INV_2ROU = pow_constexpr(TWO_ROU, FIELD_MODULUS - 2, FIELD_MODULUS);
        static const arr_u128 psi_inv_pows = get_rou_pows(INV_2ROU);
        static const bigz INV_N = pow_constexpr(N_, FIELD_MODULUS - 2, FIELD_MODULUS);
        static const bigz INV_2N = pow_constexpr(2 * N_, FIELD_MODULUS - 2, FIELD_MODULUS);
        // TODO could use intt_iter1 here? Try
        is_conv ? intt_iter(a, psi_inv_pows, INV_2N) : intt_iter1(a, psi_inv_pows, INV_N);
        a.isNTT = false;
        if (is_conv)
            assert(n == 2 * N_ && a.get(n - 1) == 0);
        return a;
    }

    using array1d<bigz, poly>::operator*;
    poly operator*(const poly& other) const {
        ASSERT(n_coeffs() == other.n_coeffs());
        return nega_ntt(other);
    }

    poly mod_cyclo(size_t N) const {
        TIMING(timing.calls_conv_to_nega += 1;)
        assert(isNTT == false);
        ASSERT(n_coeffs() == 2 * N);
        assert(get(2 * N - 1) == 0);
        poly negacyclic{N};
        for (size_t i = 0; i < N; i++)
            negacyclic.set(i, mod_sub(get_coeff(i), get_coeff(i + N)));
        return negacyclic;
    }
    G1 msm(std::vector<G1>& eval_pows_g) {
        TIMING(timing.calls_msm += 1;)
        assert(size() == 2 * N_ || size() == N_);

        // convert bigz poly coeffs to Fr scalars
        TIMING(auto start = std::chrono::high_resolution_clock::now();)
        static std::vector<Fr> scalars(2 * N_);
        #pragma omp parallel for schedule(static) num_threads(N_THREADS)
        for (size_t i = 0; i < size(); i++) {
            mpz_to_Fr(scalars[i], get(i));
        }
        TIMING(auto end = std::chrono::high_resolution_clock::now();)
        TIMING(timing.msm += end - start;)

        TIMING(auto start1 = std::chrono::high_resolution_clock::now();)
        static G1 res;
        G1::mulVecMT(res, eval_pows_g.data(), scalars.data(), size(), N_THREADS);
        TIMING(auto end1 = std::chrono::high_resolution_clock::now();)

        TIMING(timing.msm1 += end1 - start1;)
        return res;
    }
};


class hashed_t_poly : public array1d<bigz, hashed_t_poly> {
    using grp_poly = std::vector<G1>;
private:
    grp_poly arr_g; // TODO need to overload array1d destructor?
public:
    hashed_t_poly() : array1d<bigz, hashed_t_poly>() {}
    hashed_t_poly(size_t n_hashed_t_coeffs) : array1d<bigz, hashed_t_poly>(n_hashed_t_coeffs), arr_g(n_hashed_t_coeffs) {
        #pragma omp parallel for schedule(static) num_threads(N_THREADS)
        for (size_t i = 0; i < n_hashed_t_coeffs; i++) {
            set(i, 0);
            G1 x;
            x.clear();
            arr_g[i] = x;
        }
    }
    auto& get_hashed_t_coeff(size_t n) const {
        return get(n);
    }
    size_t n_hashed_t_coeffs() const {
        return size();
    }
    void set_g(size_t i, G1&& val) {
        arr_g[i] = val;
    }
    // This gets called when setting up keys. Raise (mult) Gen to get(i)
    hashed_t_poly pow() {
        hashed_t_poly res(size());
        #pragma omp parallel for schedule(static) num_threads(N_THREADS)
        for (size_t i = 0; i < size(); i++) {
            res.set_g(i, pow_(gen1, get(i)));
        }
        return res;
    }
};

// TODO eval_pows can be a poly of size 2 * N_
inline auto poly::get_hash_t(const vector_bigz& eval_pows) const {
    auto hash = get_hash(eval_pows);
    auto N = size();
    hashed_t_poly ha_poly(N);
    // TODO optimise
    vector_bigz hashed_t_vec = scalar_vec_mult(hash, eval_pows, FIELD_MODULUS);
    for (size_t i = 0; i < N; i++) {
        auto val = hashed_t_vec.at(i);
        ha_poly.set(i, val);
    }
    return ha_poly;
}

class hashed_rlwe_decomp : public array1d<bigz, hashed_rlwe_decomp> {
private:
public:
    // NOTE hashed_rlwe always has two hashed_polys
    hashed_rlwe_decomp() : array1d<bigz, hashed_rlwe_decomp>() {}
    hashed_rlwe_decomp(size_t n_hashed_polys) : array1d<bigz, hashed_rlwe_decomp>(n_hashed_polys) {}
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
    std::vector<G1> g1;
    rlwe_decomp() : array1d<poly, rlwe_decomp>() {}
    rlwe_decomp(size_t n_polys) : array1d<poly, rlwe_decomp>(n_polys), g1{n_polys} {}
    rlwe_decomp(
        size_t n_polys, size_t n_coeffs
    ) : array1d<poly, rlwe_decomp>(n_polys), g1{n_polys} {
        for (size_t i = 0; i < n_polys; i++) {
            set(i, poly(n_coeffs));
        }
    }
    poly& get_poly(size_t n) const {
        return get(n);
    }
    // calls poly's get_hash and returns hashed_rlwe_decomp
    auto get_hash(const vector_bigz& eval_pows) const {
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
    void to_eval_form() {
        for (size_t i = 0; i < n_polys(); i++) {
            set(i, get_poly(i).to_eval_form());
        }
    }
    void msm(std::vector<G1>& eval_pows_g) {
        for (size_t i = 0; i < size(); i++)
            g1.at(i) = get(i).msm(eval_pows_g);
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
    auto get_hash(const vector_bigz& eval_pows) const {
        hashed_rlwe_decomp_vec hash(n_rlwe_decomps()); // get(0).size()); // FIXME n_polys, do better (global?)
        for (size_t i = 0; i < size(); i++) {
            hash.set(i, get(i).get_hash(eval_pows));
        }
        return hash;
    }
    rlwe_decomp_vec to_eval_form() {
        for (size_t i = 0; i < n_rlwe_decomps(); i++)
            get_rlwe_decomp(i).to_eval_form();
        return *this;
    }
    void msm(std::vector<G1>& eval_pows_g) {
        for (size_t i = 0; i < size(); i++)
            get(i).msm(eval_pows_g);
    }
};

class hashed_t_rlwe : public array1d<hashed_t_poly, hashed_t_rlwe> {
private:
public:
    // NOTE hashed_t_rlwe always has two hashed_t_polys
    hashed_t_rlwe() : array1d<hashed_t_poly, hashed_t_rlwe>() {}
    hashed_t_rlwe(size_t n_hashed_t_polys) : array1d<hashed_t_poly, hashed_t_rlwe>(n_hashed_t_polys) {}
    hashed_t_rlwe(size_t n_hashed_t_polys, size_t n_coeffs) : array1d<hashed_t_poly, hashed_t_rlwe>(n_hashed_t_polys) {
        ASSERT(n_hashed_t_polys == N_POLYS_IN_RLWE);
        for (size_t i = 0; i < n_hashed_t_polys; i++) {
            set(i, hashed_t_poly(n_coeffs));
        }
    }
    auto& get_hashed_t_poly(size_t n) const {
        return get(n);
    }
    size_t n_hashed_t_polys() const {
        return size();
    }
};

class hashed_rlwe : public array1d<bigz, hashed_rlwe> {
    using int_rlwe = array1d<bigz, hashed_rlwe>;
    using grp_rlwe = std::vector<G1>;
private:
    grp_rlwe arr_g;
public:
    // NOTE hashed_rlwe always has two hashed_polys
    hashed_rlwe() : array1d<bigz, hashed_rlwe>() {}
    hashed_rlwe(size_t n_hashed_polys) : array1d<bigz, hashed_rlwe>(n_hashed_polys), arr_g(n_hashed_polys) {
        ASSERT(n_hashed_polys == N_POLYS_IN_RLWE);
        for (size_t i = 0; i < n_hashed_polys; i++) {
            G1 x;
            x.clear();
            arr_g[i] = x;
        }
    }
    size_t n_hashed_polys() const {
        return size();
    }
};

class rlwe : public array1d<poly, rlwe> {
private:
public:
    std::vector<G1> g1;
    // NOTE rlwe always has two polys
    rlwe() : array1d<poly, rlwe>() {}
    rlwe(size_t n_polys) : array1d<poly, rlwe>(n_polys), g1{n_polys} {
        ASSERT(n_polys == 2);
    }
    rlwe(size_t n_polys, size_t n_coeffs) : array1d<poly, rlwe>(n_polys), g1{n_polys} {
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
    auto get_hash(const vector_bigz& eval_pows) const {
        hashed_rlwe hash(N_POLYS_IN_RLWE);
        for (size_t i = 0; i < size(); i++) {
            hash.set(i, get(i).get_hash(eval_pows));
        }
        return hash;
    }
    // calls poly's get_hash and returns hashed_rlwe
    auto get_hash_t(const vector_bigz& eval_pows) const {
        hashed_t_rlwe hash_t(N_POLYS_IN_RLWE);
        for (size_t i = 0; i < size(); i++) {
            hash_t.set(i, get(i).get_hash_t(eval_pows));
        }
        return hash_t;
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
    rlwe_decomp decompose(const bigz& v, const u32& d, const u32& power) const {
        rlwe_decomp polys(2 * d, n_coeffs());
        for (size_t k = 0; k < N_POLYS_IN_RLWE; k++) {
            const poly& pol = get_poly(k);
            for (size_t i = 0; i < n_coeffs(); i++) {
                for (size_t j = 0; j < (size_t)d; j++) {
                    bigz decomped_coeff = (v - 1) & (pol.get_coeff(i) >> (power * j));
                    // std::cout << "decomped_coeff: " << print_to_string_mpz(decomped_coeff) << "\n";
                    ASSERT(decomped_coeff < v);
                    // polys.get_poly(d * k + j).set(i, std::move(decomped_coeff)); // XXX move doesn't make a difference here
                    polys.get_poly(d * k + j).set(i, decomped_coeff);
                }
            }
        }
        return polys;
    }
    auto mod_cyclo(size_t N) const {
        rlwe conv(N_POLYS_IN_RLWE);
        for (size_t i = 0; i < N_POLYS_IN_RLWE; i++) {
            poly p = get_poly(i).mod_cyclo(N);
            conv.set(i, p);
        }
        return conv;
    }
    void to_eval_form() {
        for (size_t i = 0; i < N_POLYS_IN_RLWE; i++) {
            set(i, get_poly(i).to_eval_form());
        }
    }
    void to_coeff_form() {
        for (size_t i = 0; i < N_POLYS_IN_RLWE; i++) {
            set(i, get_poly(i).to_coeff_form());
        }
    }
    void msm(std::vector<G1>& eval_pows_g) {
        for (size_t i = 0; i < size(); i++) {
            g1.at(i) = get(i).msm(eval_pows_g);
        }
        // std::cout << "finished 1 rlwe\n";
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
    auto get_hash(const vector_bigz& eval_pows) const {
        hashed_rlwe_vec hash(size());
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
    auto decompose(const bigz& v, u32 d, u32 power) const {
        rlwe_decomp_vec decomps(size());//, 2 * d, n_coeffs());
        for (size_t i = 0; i < size(); i++)
            decomps.set(i, get_rlwe(i).decompose(v, d, power));
        return decomps;
    }
    auto mod_cyclo(size_t N) const {
        rlwe_vec conv(n_rlwes());
        for (size_t i = 0; i < n_rlwes(); i++) {
            rlwe r = get_rlwe(i).mod_cyclo(N);
            conv.set(i, r);
        }
        return conv;
    }
    void to_coeff_form() const {
        for (size_t i = 0; i < n_rlwes(); i++) {
            get_rlwe(i).to_coeff_form();
        }
    }
    void msm(std::vector<G1>& eval_pows_g) {
        for (size_t i = 0; i < size(); i++) {
            get(i).msm(eval_pows_g);
        }
    }
};

class hashed_t_rgsw : public array1d<bigz, hashed_t_rgsw> {
private:
std::vector<G2> g2;
public:
    using hashed_t_poly = bigz;
    hashed_t_rgsw() : array1d<hashed_t_poly, hashed_t_rgsw>() {}
    hashed_t_rgsw(size_t n_hashed_t_polys) : array1d<hashed_t_poly, hashed_t_rgsw>(n_hashed_t_polys), g2{n_hashed_t_polys} {
    }
    auto& get_hashed_t_poly(size_t n) const {
        return get(n);
    }
    size_t n_hashed_t_polys() const {
        return size();
    }
    GT get_hash_sec(const rlwe_decomp& other) const {
        assert(size() == other.size());
        GT res{1};
        for (size_t i = 0; i < n_hashed_t_polys(); i++) {
            GT tmp;
            pairing(tmp, other.g1[i], g2[i]);
            res *= tmp;
        }
        return res;
    }
    hashed_t_rgsw pow() const {
        hashed_t_rgsw x(n_hashed_t_polys());
        for (size_t i = 0; i < n_hashed_t_polys(); i++) {
            x.g2[i] = pow2_(gen2, get(i));
        }
        return x;
    }
};

class hashed_rgsw : public array1d<bigz, hashed_rgsw> {
    using grp_rgsw = std::vector<G1>;
private:
    grp_rgsw arr_g;
public:
    hashed_rgsw() : array1d<bigz, hashed_rgsw>() {}
    hashed_rgsw(
        size_t n_hashed_polys
    ) : array1d<bigz, hashed_rgsw>(n_hashed_polys), arr_g(n_hashed_polys) {
        for (size_t i = 0; i < n_hashed_polys; i++) {
            set(i, 0);
            G1 x;
            x.clear();
            arr_g[i] = x;
        }
    }
    G1 get_g(size_t i) const {
        return arr_g[i];
    }
    void set_g(size_t i, G1 val) {
        arr_g[i] = val;
    }

    auto& get_hashed_poly(size_t n) const {
        return get(n);
    }
    size_t n_hashed_polys() const {
        return size();
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

    rlwe convolve(const rlwe_decomp& other) const {
        ASSERT(n_rlwes() == other.n_polys());
        rlwe res(n_polys(), 2 * N_); // XXX 2 * N_ zero-inits polys for += later
        for (size_t i = 0; i < n_rlwes(); i++) {
            for (size_t j = 0; j < n_polys(); j++) {
                res.get(j) += get(i).get(j).convolve(other.get(i));
            }
        }
        return res;
    }

    void to_eval_form() {
        for (size_t i = 0 ; i < n_rlwes(); i++)
            get_rlwe(i).to_eval_form();
    }
};


class flat_rgsw : public array1d<poly, flat_rgsw> {
private:
public:
    flat_rgsw() : array1d<poly, flat_rgsw>() {}
    flat_rgsw(size_t n_polys) : array1d<poly, flat_rgsw>(n_polys) {}
    flat_rgsw(
        size_t n_polys, size_t n_coeffs
    ) : array1d<poly, flat_rgsw>(n_polys) {
        for (size_t i = 0; i < n_polys; i++)
            set(i, poly(n_coeffs));
    }
    poly& get_poly(size_t n) const {
        return get(n);
    }
    void set_poly(size_t n, const poly& val) {
        set(n, val);
    }
    size_t n_polys() const {
        return size();
    }
    size_t n_coeffs() const {
        return get_poly(0).n_coeffs();
    }

    // calls rlwe's get_hash for all the rlwe's in (*this)
    auto get_hash(const vector_bigz& eval_pows) const {
        hashed_rgsw hash(n_polys()); // FIXME n_polys, do better (global?)
        for (size_t i = 0; i < n_polys(); i++) {
            hash.set(i, get(i).get_hash(eval_pows));
        }
        return hash;
    }
    // calls rlwe's get_hash_t for all the rlwe's in (*this)
    auto get_hash_t(const vector_bigz& eval_pows) const {
        hashed_t_rgsw hash(n_polys()); // , n_coeffs()); // FIXME n_polys, do better (global?)
        for (size_t i = 0; i < size(); i++) {
            hash.set(i, get(i).get_hash(eval_pows));
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
        size_t n_hashed_rgsws, size_t n_hashed_polys
    ) : array1d<hashed_rgsw, hashed_rgsw_vec>(n_hashed_rgsws) {
        for (size_t i = 0; i < n_hashed_rgsws; i++)
            set(i, hashed_rgsw(n_hashed_polys));
    }
    hashed_rgsw& get_hashed_rgsw(size_t n) const {
        return get(n);
    }
    size_t n_hashed_rgsws() const {
        return size();
    }
    size_t n_hashed_polys() const {
        size_t n_hashed_polys = get_hashed_rgsw(0).n_hashed_polys();
        return n_hashed_polys;
    }
    bigz dot_prod(const hashed_rlwe_decomp_vec& rlwe_v) const {
        ASSERT(size() == rlwe_v.size());
        ASSERT(get(0).size() == rlwe_v.get(0).size());
        bigz sum = 0;
        for (size_t i = 0; i < n_hashed_rgsws(); i++) {
            for (size_t j = 0; j < n_hashed_polys(); j++) {
                sum += get(i).get(j) * rlwe_v.get(i).get(j);
            }
        }
        return mod_(sum, FIELD_MODULUS);
    };
};

class hashed_t_rgsw_vec : public array1d<hashed_t_rgsw, hashed_t_rgsw_vec> {
private:
public:
    hashed_t_rgsw_vec() : array1d<hashed_t_rgsw, hashed_t_rgsw_vec>() {}
    hashed_t_rgsw_vec(size_t n_hashed_t_rgsws) : array1d<hashed_t_rgsw, hashed_t_rgsw_vec>(n_hashed_t_rgsws) {}
    hashed_t_rgsw& get_hashed_t_rgsw(size_t n) const {
        return get(n);
    }
    size_t n_hashed_t_rgsws() const {
        return size();
    }
    size_t n_hashed_t_polys() const {
        return get_hashed_t_rgsw(0).size();
    }
    GT get_hash_sec(const rlwe_decomp_vec& other) const {
        ASSERT(n_hashed_t_rgsws() == other.n_rlwe_decomps());
        GT res{1};
        for (size_t i = 0; i < n_hashed_t_rgsws(); i++)
            res *= get(i).get_hash_sec(other.get(i));
        return res;
    }
    hashed_t_rgsw_vec pow() const {
        hashed_t_rgsw_vec x(n_hashed_t_rgsws());
        for (size_t i = 0; i < n_hashed_t_rgsws(); i++) {
            x.set(i, get_hashed_t_rgsw(i).pow());
        }
        return x;
    }
};


class flat_rgsw_vec : public array1d<flat_rgsw, flat_rgsw_vec> {
private:
public:
    flat_rgsw_vec(size_t n_flat_rgsws) : array1d<flat_rgsw, flat_rgsw_vec>(n_flat_rgsws) {}
    flat_rgsw_vec(
        size_t n_flat_rgsws, size_t n_polys, size_t n_coeffs
    ) : array1d<flat_rgsw, flat_rgsw_vec>(n_flat_rgsws) {
        for (size_t i = 0; i < n_flat_rgsws; i++)
            set(i, flat_rgsw(n_polys, n_coeffs));
    }
    flat_rgsw& get_flat_rgsw(size_t n) const {
        return get(n);
    }
    size_t n_flat_rgsws() const {
        return size();
    }
    size_t n_polys() const {
        return get_flat_rgsw(0).n_polys();
    }
    size_t n_coeffs() const {
        return get_flat_rgsw(0).n_coeffs();
    }
    // calls flat_rgsw's get_hash
    auto get_hash(const vector_bigz& eval_pows) const {
        hashed_rgsw_vec hash(n_flat_rgsws());
        for (size_t i = 0; i < n_flat_rgsws(); i++) {
            hash.set(i, get(i).get_hash(eval_pows));
        }
        return hash;
    }
    // calls flat_rgsw's get_hash_t
    auto get_hash_t(const vector_bigz& eval_pows) const {
        hashed_t_rgsw_vec hash(n_flat_rgsws());
        for (size_t i = 0; i < size(); i++) {
            hash.set(i, get(i).get_hash_t(eval_pows));
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

    rlwe_vec convolve(const rlwe_decomp_vec& other) const {
        ASSERT(n_cols() == other.size());
        ASSERT(n_coeffs() == 2 * N_);
        rlwe_vec res(n_rows());
        for (size_t i = 0; i < n_rows(); i++) {
            rlwe sum(n_polys(), 2 * N_); // FIXME
            for (size_t j = 0; j < n_cols(); j++) {
                sum += get(i, j).convolve(other.get(j));
            }
            res.set(i, sum);
        }
        return res;
    }
    rlwe_vec convolve(rlwe_decomp_vec&& other) const {
        return convolve(other);
    }
    void to_eval_form() {
        for (size_t row = 0; row < n_rows(); row++) {
            for (size_t col = 0; col < n_cols(); col++) {
                get_rgsw(row, col).to_eval_form();
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

class hashed_t_veri_vec_inner : public array1d<bigz, hashed_t_veri_vec_inner> {
private:
    G2 g2[2];
public:
    using htvvi = hashed_t_veri_vec_inner;
    using htvvi_arr = array1d<bigz, hashed_t_veri_vec_inner>;
    hashed_t_veri_vec_inner() : htvvi_arr() {}
    hashed_t_veri_vec_inner(size_t n_ha_polys) : htvvi_arr{n_ha_polys} {
        assert(n_ha_polys == N_POLYS_IN_RLWE);
    }

    // XXX pairing
    GT get_hash_sec(const rlwe& other) const {
        assert(size() == other.size());
        GT res{1};
        for (size_t i = 0; i < size(); i++) {
            GT tmp;
            pairing(tmp, other.g1.at(i), g2[i]);
            res *= tmp;
        }
        return res;
    }

    // XXX
    htvvi pow() const {
        htvvi res{size()};
        for (size_t i = 0; i < size(); i++) {
            res.g2[i] = pow2_(gen2, get(i));
        }
        return res;
    }
};

class hashed_t_veri_vec : public array1d<hashed_t_veri_vec_inner, hashed_t_veri_vec> {
public:
    using htvv = hashed_t_veri_vec;
    using htvv_arr = array1d<hashed_t_veri_vec_inner, hashed_t_veri_vec>;
    hashed_t_veri_vec() : htvv_arr() {}
    hashed_t_veri_vec(size_t n_inner) : htvv_arr(n_inner) {}

    GT get_hash_sec(const rlwe_vec& other) const {
        assert(size() == other.size());
        GT res{1};
        for (size_t i = 0; i < size(); i++) {
            res *= get(i).get_hash_sec(other.get(i));
        }
        return res;
    }
    // TODO pow can be made a memb func on array1d again - recurse to ha_poly
    htvv pow() const {
        htvv res{size()};
        for (size_t i = 0; i < size(); i++) {
            res.set(i, get(i).pow());
        }
        return res;
    }
};

class veri_vec_inner : public array1d<bigz, veri_vec_inner> {
public:
using vvi_arr = array1d<bigz, veri_vec_inner>;
    bigz dot_prod(const hashed_rlwe& other) const {
        assert(size() == other.size());
        bigz res{0};
        for (size_t i = 0; i < size(); i++) {
            bigz val = other.get(i) * get(i);
            val = mod_(val, FIELD_MODULUS);
            res += val;
        }
        return res;
    }
    using vvi_arr::operator*;
    flat_rgsw operator*(const rgsw& other) {
        assert(size() == other.n_polys());
        flat_rgsw res{other.n_rlwes()};
        for (size_t i = 0; i < other.size(); i++) { // columns (rlwe's of rgsw)
            poly a{N_}; // TODO can grab el 0 then += from idx 1 onwards
            for (size_t j = 0; j < size(); j++) { // 2 bigz (vvi), 2 polys (rgsw)
                a += other.get(i).get(j) * get(j);
            }
            res.set(i, a);
        }
        return res;
    }
    using htvvi = hashed_t_veri_vec_inner;

    // TODO remove, should be copy
    htvvi get_hash_t() const {
        htvvi res{N_POLYS_IN_RLWE};
        for (size_t i = 0; i < size(); i++) {
            res[i] = get(i);
        }
        return res;
    }
};

class veri_vec : public array1d<veri_vec_inner, veri_vec> {
public:
using vv_arr = array1d<veri_vec_inner, veri_vec>;
    using vv_arr::operator*;
    flat_rgsw_vec operator*(const rgsw_mat& other) const {
        assert(size() == other.n_rows());
        flat_rgsw_vec res(other.n_cols());
        for (size_t j = 0; j < other.n_cols(); j++) {
            // flat_rgsw el{get(0) * other.get(0, j)}; // FIXME copy ctor issue
            flat_rgsw el = get(0) * other.get(0, j);
            for (size_t i = 1; i < size(); i++) {
                el += get(i) * other.get(i, j);
            }
            res.set(j, el);
        }
        return res;
    }
    bigz dot_prod(const hashed_rlwe_vec& other) const {
        // does dot prod, calling inner's dot prod
        assert(size() == other.size());
        bigz res{0}; // TODO can grab el 0 then += from idx 1 onwards
        for (size_t i = 0; i < size(); i++)
            res += get(i).dot_prod(other.get(i));
        return res;
    }
    using htvv = hashed_t_veri_vec;
    // TODO remove, should be copy
    htvv get_hash_t() const {
        htvv res{size()}; // FIXME use {} or ()?
        // std::cout << "htvv::size(): " << res.size() << "\n";
        for (size_t i = 0; i < size(); i++) {
            res.set(i, get(i).get_hash_t());
        }
        return res;
    }
};

bigz random_i128(bigz from, bigz to_inclusive) {
    // random int in range [from, to_inclusive]
    // TODO docs suggest the above inclusive range. Double check
    // HACK
    long int from_ = mpz_get_si(from.get_mpz_t());
    long int to_inclusive_ = mpz_get_si(to_inclusive.get_mpz_t());
    std::random_device rd;  // Non-deterministic seed
    std::mt19937 gen(rd()); // Mersenne Twister engine
    std::uniform_int_distribution<unsigned long> distrib(from_, to_inclusive_); // Range: 0 to 100
    bigz random_number = distrib(gen);
    return random_number;
}

u32 get_decomp_power() {
    // v - power of 2, s.t. v^{d-1} < q < v^d
    double log2mpz = log2_mpz(FIELD_MODULUS);
    size_t q_bits = static_cast<size_t>(log2mpz) + 1;
    u32 remainder = q_bits % N_DECOMP;
    u32 power = remainder ? q_bits / N_DECOMP + 1 : q_bits / N_DECOMP;
    assert(FIELD_MODULUS <= bigz(1) << (power * N_DECOMP));
    assert(FIELD_MODULUS > bigz(1) << (power * (N_DECOMP - 1)));
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

    // ======== Scale up G, R, and H to integers ========
    array2d<bigz> scalar_mat_mult(double scalar, matrix_double& mat, bigz q, size_t rows, size_t cols) {
        array2d<bigz> result(rows, cols);
        for (size_t i = 0; i < rows; i++) {
            for (size_t j = 0; j < cols; j++) {
                mpf_class rounded = mpf_round(scalar * mat[i][j]);
                bigz modded = (bigz)rounded;
                modded = mod_(modded, q);
                result.set(i, j, modded);
            }
        }
        return result;
    }

    // FIXME extract lambdas from vfhe.cpp and remove this double up
    // TODO update to array1d<bigz> (add as another class?)
    vector_bigz scalar_vec_mult(double scalar, vector_double& vec, bigz q) {
        // TODO update to array1d<bigz> (add as another class?)
        vector_bigz result(vec.size());
        ASSERT(result.capacity() == result.size());
        for (size_t i = 0; i < vec.size(); i++) {
            mpf_class rounded = mpf_round(scalar * vec[i]);
            bigz modded = (bigz)rounded;
            modded = mod_(modded, q);
            result.at(i) = modded;
        }
        return result;
    }

public:
    Params() :
        N(N_),
        iter_(3),
        s(10000.0),
        L(10000.0),
        r(10000.0),
        q(FIELD_MODULUS),
        d{N_DECOMP},
        power{get_decomp_power()},
        v{static_cast<bigz>(1) << power},
        F(F_.size(), F_.at(0).size()),
        G_bar(scalar_mat_mult(s, G, q, G.size(), G.at(0).size())),
        R_bar(scalar_mat_mult(s, R, q, R.size(), R.at(0).size())),
        H_bar(scalar_mat_mult(s, H, q, H.size(), H.at(0).size())),
        x_plant_init{1, -1, 0, 0.7, 2},
        x_cont_init{-0.001, 0.013, 0.2, -0.02, 0},
        x_cont_init_scaled{scalar_vec_mult(r * s * L, x_cont_init, q)}
    {
        // TODO update or remove
        ASSERT(A.size() == B.size());
        ASSERT(A.size() == F_.size());
        assert(d >= 1);
        // Construct F from F_
        for (size_t i = 0; i < F_.size(); i++) {
            for (size_t j = 0; j < F_.at(0).size(); j++) {
                mpf_class val = F_.at(i).at(j);
                bigz val1 = bigz(val); // TODO - mpf_class intermediate required?
                if (val1 < 0)
                    val1 + FIELD_MODULUS;
                F.set(i, j, val1);
            }
        }
    }
    size_t N, iter_;
    double s, L, r;
    bigz q;
    u32 d;
    u32 power;
    bigz v;
    array2d<bigz> F, G_bar, R_bar, H_bar;

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
    vector_double x_plant_init;
    vector_double x_cont_init;
    vector_bigz x_cont_init_scaled;
    void print() {
        std::cout << "*** START Params.print ***\n";
        std::cout << "N_THREADS: " << N_THREADS << "\n";
        std::cout << "N: " << N << "\n";
        std::cout << "q: " << print_to_string_mpz(q) << "\n";
        std::cout << "d: " << d << "\n";
        std::cout << "v: " << v << "\n";
        std::cout << "power: " << power << "\n";
        std::cout << "s: " << s << "\n";
        std::cout << "L: " << L << "\n";
        std::cout << "r: " << r << "\n";
        std::cout << "iter: " << iter_ << "\n";
        // std::cout << "F:\n";
        // F.print();
        // std::cout << "G_bar:\n";
        // G_bar.print();
        // std::cout << "R_bar:\n";
        // R_bar.print();
        // std::cout << "H_bar:\n";
        // H_bar.print();
        // std::cout << "x_cont_init_scaled:\n";
        // print_vector_mpz(x_cont_init_scaled);
        std::cout << "*** END Params.print ***\n\n\n";
    }
    // TODO pass in gen as a param
    vector_bigz sample_knowledge_exponents(bigz from, bigz to_inclusive) {
        size_t N = 6;
        vector_bigz res(N);
        for (size_t i = 0; i < N; i++) {
            // TODO update range
            res.at(i) = random_i128(from, to_inclusive);
            DEBUG1(res.at(i) = 1;)
        }
        // res <- {alpha_0, alpha_1, gamma_0, gamma_1, rho_0, rho_1}
        return res;
    }

    // TODO pass in gen as a param
    std::vector<veri_vec> sample_verification_vectors(__uint128_t m, __uint128_t n, bigz from, bigz to_inclusive) {
        size_t N = 3;
        std::vector<veri_vec> res{N};
        for (size_t i = 0; i < N - 1; i++) {
            res.at(i) = veri_vec{n};
        }
        res.at(N - 1) = veri_vec{m};

        for (auto& x : res) {
            for (size_t i = 0; i < x.size(); i++) {
                x.set(i, veri_vec_inner{N_POLYS_IN_RLWE});
                for (size_t j = 0; j < N_POLYS_IN_RLWE; j++) {
                    x.get(i).set(j, random_i128(from, to_inclusive));
                    DEBUG1(x.at(i) = 1;)
                }
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
    hashed_t_veri_vec gr;
    hashed_t_rgsw_vec grFr;
    hashed_t_rgsw_vec gsHr;
    hashed_t_veri_vec gr_rho;
    hashed_t_rgsw_vec grFr_alpha;
    hashed_t_rgsw_vec gsHr_gamma;
};

struct veri_key {
    veri_vec s;
    hashed_rgsw_vec rG_0;
    hashed_rgsw_vec rG_1;
    hashed_rgsw_vec rR_0;
    hashed_rgsw_vec rR_1;
    bigz rho_0;
    bigz rho_1;
    bigz alpha_0;
    bigz alpha_1;
    bigz gamma_0;
    bigz gamma_1;
};

struct check_key {
    flat_rgsw_vec r_1_rgsw;
    flat_rgsw_vec r_0_rgsw;
};

struct Proof {
    // GT g_1;
    GT grx_;
    GT grx_hi;
    GT grFrx;
    GT gsHrx;
    GT gr_rho_x_;
    GT gr_rho_x_hi;
    GT grFr_alpha_x;
    GT gsHr_gamma_x;
};