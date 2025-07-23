// TODO write more tests
// TODO rgsw code
    // TODO PRNG
    // TODO decomposition
    // TODO enc/dec
// TODO control code
// TODO verification code
// TODO proof code (linear/dynamic checks)
// TODO NTT code
// TODO homomorphic hash variant
// TODO parallellise
// TODO cyclic group code
// TODO replace *this and (*this).member_func() with this and this->member_func() ???
/*
Setup:
rF, rG, sH
- veri_vec * rgsw_vec
    - scalar * rgsw

Verfication:
rx' = rFx + rGx
su = sHu
*/

#include <iostream>
#include <cstdlib>
#include <vector>
#include <cassert>
#include <cmath>
#include <climits>
#include <algorithm>
#include <type_traits>
#include <chrono>
#include <random>
#include "omp.h"

typedef __uint128_t i128;
// constexpr i128 GROUP_MODULUS = 23;
// constexpr i128 FIELD_MODULUS = 11;
// constexpr i128 GENERATOR = 4;
constexpr i128 GROUP_MODULUS = 540431955285196831;
constexpr i128 FIELD_MODULUS = 18014398509506561;
constexpr i128 GENERATOR = 1073741824;

// FIXME make types __uint128 so that regular modding works
i128 mod_(__int128_t val, i128 q) {
    val %= q;
    if (val < 0) {
        val = (val + q) % q;
    }
    return val;
}

std::vector<i128> sample_discrete_gaussian(size_t N, double mu = 3.2, double sigma = 19.2) {
    std::vector<i128> result(N);
    std::random_device rd;
    std::mt19937 gen(rd());
    std::normal_distribution<double> dist(mu, sigma);
    for (size_t i = 0; i < N; ++i) {
        result[i] = static_cast<i128>(std::round(dist(gen)));
    }
    return result;
}
std::vector<i128> sample_secret_key(size_t N) {
    // Sample a secret key for the RGSW scheme.
    // Each entry is -1, 0, or 1, with probabilities 0.25, 0.5, 0.25 respectively.
    std::vector<i128> s(N);
    std::random_device rd;
    std::mt19937 gen(rd());
    std::discrete_distribution<int> dist({0.25, 0.5, 0.25});
    for (size_t i = 0; i < N; ++i) {
        int val = dist(gen);
        if (val == 0) s[i] = -1;
        else if (val == 1) s[i] = 0;
        else s[i] = 1;
    }
    return s;
}
// Sample a random polynomial of degree N-1 with coefficients in the range [0, q).
std::vector<i128> sample_random_polynomial(size_t N, i128 q) {
    std::vector<i128> poly(N);
    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_int_distribution<i128> dist(0, q - 1);
    for (size_t i = 0; i < N; ++i) {
        poly[i] = dist(gen);
    }
    return poly;
}

// Sample a noise polynomial of degree N-1 with coefficients from the discrete Gaussian distribution.
std::vector<i128> sample_noise_polynomial(size_t N, double mu = 3.2, double sigma = 19.2) {
    return sample_discrete_gaussian(N, mu, sigma);
}

std::string print_to_string_i128(i128 n) {
    if (n == 0) {
        return "0";
    }
    bool neg = false;
    // if (n < 0) {
    //     neg = true;
    //     n = -n;
    // }
    std::string buf;
    while (n > 0) {
        buf += '0' + (n % 10);
        n /= 10;
    }
    if (neg) buf += '-';
    std::reverse(buf.begin(), buf.end());
    return buf;
}

template<typename T, typename Derived>
class array1d {
private:
    size_t size_;
    T* arr;
public:
    array1d() : size_(0) {}
    array1d(size_t size) : size_(size) {
        arr = new T[size]();
    }
    ~array1d() = default;
    void check_index_bounds(size_t n) const {
        if (n >= size_) {
            throw std::out_of_range(
                "Index error: accessing arr[" + std::to_string(n) + "]"
                + " in a " + std::to_string(size_) + " element array."
            );
        }
    }
    T& get(size_t n) const {
        // check_index_bounds(n);
        return arr[n];
    }
    void check_value_bounds(T val) {
        if constexpr (std::is_same_v<T, i128>) {
            i128 min_val = -1;
            min_val <<= 63;
            i128 max_val = (1UL << 63) - 1;
            if (val < min_val || val > max_val) {
                throw std::out_of_range(
                    "(array1d) Value out of range: " + print_to_string_i128(val)
                );
            }
        }
    }
    void set(int n, T val) {
        // check_index_bounds(n);
        // check_value_bounds(val);
        arr[n] = val;
    }
    size_t size() const {
        return size_;
    }
    Derived operator*(const Derived& other) const {
        size_t N = size();
        Derived res(N);
        for (size_t i = 0; i < N; i++) {
            auto val = (get(i) * other.get(i));
            if constexpr (std::is_same_v<T, i128>)
                val %= FIELD_MODULUS;
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
                val %= FIELD_MODULUS;
            res.set(i, val);
        }
        return res;
    }
    Derived group_mult(const Derived& other) const {
        size_t N = size();
        Derived res(N);
        for (size_t i = 0; i < N; i++) {
            T val;
            if constexpr (std::is_same_v<T, i128>) {
                val = (get(i) * other.get(i));
                val %= GROUP_MODULUS;
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
                val = get(i) * scalar;
                val %= GROUP_MODULUS;
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
                val %= FIELD_MODULUS;
            res.set(i, val);
        }
        return res;
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
        power %= FIELD_MODULUS;
        base %= GROUP_MODULUS;
        i128 result = 1;
        i128 one = 1;
        i128 two = 2;
        while (power > 0) {
            if ((power % two) == one)
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

    void print() {
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
    array2d() : rows_(0), cols_(0) {}
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
        // check_index_bounds(row, col);
        return arr[row][col];
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
        }
    }
    void set(size_t row, size_t col, T val) {
        // check_index_bounds(row, col);
        // check_value_bounds(val);
        arr[row][col] = val;
    }
    size_t n_rows() const {
        return rows_;
    }
    size_t n_cols() const {
        return cols_;
    }
    void print() const {
        for (size_t i = 0; i < rows_; i++) {
            std::cout << "{";
            for (size_t j = 0; j < cols_ - 1; j++) {
                std::cout << print_to_string_i128(arr[i][j]) << ", ";
            }
            std::cout << print_to_string_i128(arr[i][cols_ - 1]) << "}\n";
        }
    }
};

class poly : public array1d<i128, poly> {
private:
    bool isNTT = false;
public:
    poly() : array1d<i128, poly>() {}
    poly(size_t N) : array1d<i128, poly>(N) {
        for (size_t i = 0; i < N; i++) {
            set(i, 0);
        }
    }
    auto& get_coeff(size_t n) const {
        return get(n);
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
    void set_coeffs_to_one() {
        for (auto& p : *this) {
            for (size_t j = 0; j < p.size(); j++) {
                p.set(j, 1);
            }
        }
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
    size_t get_n_polys() const {
        return get_rlwe_decomp(0).size();
    }
    size_t get_n_coeffs() const {
        return get_rlwe_decomp(0).get_poly(0).size();
    }
};


class rgsw : public array1d<rlwe, rgsw> {
private:
public:
    using array1d<rlwe, rgsw>::operator*;
    using array1d<rlwe, rgsw>::pow;
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
    size_t get_n_polys() const {
        return get_rlwe(0).size();
    }
    size_t get_n_coeffs() const {
        return get_rlwe(0).get_poly(0).size();
    }
    rlwe operator*(const rlwe_decomp& other) const {
        size_t N = size();
        assert(N == other.size());
        rlwe res(get_n_polys(), get_n_coeffs());
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

    rlwe pow(const rlwe_decomp& other) const {
        size_t n_polys = get_n_polys();
        size_t n_coeffs = get_n_coeffs();
        size_t N = size();
        assert(N == other.size());
        rlwe_vec res_vec(N, n_polys, n_coeffs);
        #pragma omp parallel for num_threads(2)
        for (size_t i = 0; i < N; i++) {
            // int thread_id = omp_get_thread_num();
            // #pragma omp critical
            // {
            //     std::cout << "pow1: hi from thread #" << thread_id << std::endl;
            // }

            rlwe& res = res_vec.get(i);
            auto val0 = get(i).get(0).pow(other.get(i));
            auto val1 = get(i).get(1).pow(other.get(i));
            res.set(0, val0);
            res.set(1, val1);
        }

        rlwe res(n_polys, n_coeffs);
        res.set_coeffs_to_one();
        for (auto& rlwe_ : res_vec) {
            res = res.group_mult(rlwe_);
        }
        return res;
    }
};

class rgsw_vec : public array1d<rgsw, rgsw_vec> {
private:
public:
    using array1d<rgsw, rgsw_vec>::pow;
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
    size_t get_n_rlwes() const {
        return get_rgsw(0).size();
    }
    size_t get_n_polys() const {
        return get_rgsw(0).get_n_polys();
    }
    size_t get_n_coeffs() const {
        return get_rgsw(0).get_n_coeffs();
    }
    rlwe operator*(const rlwe_decomp_vec& other) const {
        size_t n = other.size();
        assert(size() == n);
        rlwe sum(get_n_polys(), get_n_coeffs());
        for (size_t i = 0; i < n; i++) {
            auto val = get(i) * other.get(i);
            sum = sum + val;
        }
        return sum;
    }
    rlwe pow(const rlwe_decomp_vec& other) const {
        size_t n = other.size();
        assert(size() == n);
        size_t n_polys = get_n_polys();
        size_t n_coeffs = get_n_coeffs();
        rlwe_vec res_vec(n, n_polys, n_coeffs);
        rlwe res(n_polys, n_coeffs);
        res.set_coeffs_to_one();
        // size_t n_threads_max = omp_get_max_threads();
        // std::cout << "Max threads: " << n_threads_max << std::endl;
        #pragma omp parallel for num_threads(12)
        for (size_t i = 0; i < n; i++) {
            // int thread_id = omp_get_thread_num();
            // #pragma omp critical
            // {
            //     std::cout << "pow0: hi from thread #" << thread_id << std::endl;
            // }
            auto val = get(i).pow(other.get(i));
            res_vec.set(i, val);
        }
        for (auto& prod : res_vec) {
            res = res.group_mult(prod);
        }
        return res;
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
    ~rgsw_mat() {}

    rlwe_vec operator*(const rlwe_decomp_vec& other) const {
        size_t rows = n_rows();
        size_t cols = n_cols();
        assert(cols == other.size());
        size_t n_polys = get_n_polys();
        size_t n_coeffs = get_n_coeffs();
        size_t n_coeffs_other = other.get_n_coeffs();
        assert(n_coeffs == n_coeffs_other);

        rlwe_vec res(rows, n_polys, n_coeffs);
        for (size_t i = 0; i < rows; i++) {
            rlwe sum(n_polys, n_coeffs);
            for (size_t j = 0; j < cols; j++) {
                auto val = get(i, j) * other.get(j);
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
        size_t n_coeffs = rg.get_n_coeffs();
        size_t n_polys = rg.get_n_polys();
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
    size_t get_n_rlwes() const {
        return get_rgsw(0, 0).size();
    }
    size_t get_n_polys() const {
        return get_rgsw(0, 0).get_n_polys();
    }
    size_t get_n_coeffs() const {
        return get_rgsw(0, 0).get_n_coeffs();
    }
};


class veri_vec_scalar : public array1d<i128, veri_vec_scalar> {
private:
public:
    using array1d<i128, veri_vec_scalar>::pow;
    veri_vec_scalar() : array1d<i128, veri_vec_scalar>() {}
    veri_vec_scalar(size_t N) : array1d<i128, veri_vec_scalar>(N) {}
    ~veri_vec_scalar() {}

    rgsw_vec operator*(const rgsw_mat& other) const {
        size_t rows = other.n_rows();
        size_t cols = other.n_cols();
        assert(rows == size());
        size_t n_rlwes = other.get_n_rlwes();
        size_t n_polys = other.get_n_polys();
        size_t n_coeffs = other.get_n_coeffs();
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

void init(rgsw_mat& F, rlwe_decomp_vec& x, veri_vec_scalar& r) {
    std::random_device rd;  // Non-deterministic seed
    std::mt19937 gen(rd()); // Mersenne Twister engine
    std::uniform_int_distribution<i128> distrib(0, FIELD_MODULUS - 1); // Range: 0 to 100
    i128 random_number = distrib(gen);
    i128 counter = random_number;
    for (size_t i = 0; i < F.n_rows(); ++i) {
        for (size_t j = 0; j < F.n_cols(); ++j) {
            rgsw& rg = F.get_rgsw(i, j);
            for (rlwe& rl : rg) {
                for (poly& p : rl) {
                    for (size_t m = 0; m < p.size(); ++m) {
                        p.set(m, counter++ % FIELD_MODULUS); // just an example initialization
                    }
                }
            }
        }
    }
    counter = distrib(gen);
    for (size_t i = 0; i < x.size(); ++i) {
        rlwe_decomp& rd = x.get(i);
        for (size_t j = 0; j < rd.size(); ++j) {
            poly& p = rd.get(j);
            for (size_t m = 0; m < p.size(); ++m) {
                p.set(m, counter++ % FIELD_MODULUS);
            }
        }
    }
    counter = distrib(gen);
    for (size_t i = 0; i < r.size(); ++i) {
        r.set(i, counter++ % FIELD_MODULUS);
    }
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
                result.set(i, j, mod_(static_cast<__int128_t>(std::round(scalar * mat[i][j])), q));
            }
        }
        return result;
    }

    // TODO update to array1d<i128> (add as another class?)
    poly scalar_vec_mult(i128 scalar, std::vector<double>& vec, i128 q) {
        // TODO update to array1d<i128> (add as another class?)
        poly result(vec.size());
        for (size_t i = 0; i < vec.size(); ++i) {
            result.set(i, mod_(static_cast<__int128_t>(std::round(scalar * vec[i])), q));
        }
        return result;
    }

    void generate_field_and_group_params() {
        // TODO use NTL and add funcs from Python implementation
        p = 540431955285196831;
        q = 18014398509506561;
        g = 1073741824;
    }

public:
    Params() :
        p(540431955285196831),
        q(18014398509506561),
        g(1073741824),
        s(10000),
        L(10000),
        r(10000),
        iter_(100),
        F(scalar_mat_mult(s, F_, q, F_.size(), F_.at(0).size())),
        G_bar(scalar_mat_mult(s, G, q, G.size(), G.at(0).size())),
        R_bar(scalar_mat_mult(s, R, q, R.size(), R.at(0).size())),
        H_bar(scalar_mat_mult(s, H, q, H.size(), H.at(0).size()))
    {
        // generate_field_and_group_params();
        x_cont_init_scaled = scalar_vec_mult(r * s * L, x_cont_init, q);
        // TODO update or remove
        assert(A.size() == B.size());
        assert(A.size() == F_.size());
    }
    i128 p, q, g, s, L, r, iter_, x_dim, y_dim, u_dim;
    array2d<i128> F, G_bar, R_bar, H_bar;
    // TODO update to array1d<i128> (add as another class?)
    poly x_cont_init_scaled;

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
        // print these: p, q, g, s, L, r, iter_;
        std::cout << "p: " << print_to_string_i128(p) << "\n";
        std::cout << "q: " << print_to_string_i128(q) << "\n";
        std::cout << "g: " << print_to_string_i128(g) << "\n";
        std::cout << "s: " << print_to_string_i128(s) << "\n";
        std::cout << "L: " << print_to_string_i128(L) << "\n";
        std::cout << "r: " << print_to_string_i128(r) << "\n";
        std::cout << "iter: " << print_to_string_i128(iter_) << "\n";
        std::cout << "G_bar:\n";
        G_bar.print();
        std::cout << "R_bar:\n";
        R_bar.print();
        std::cout << "H_bar:\n";
        H_bar.print();
        std::cout << "x_cont_init_scaled:\n";
        x_cont_init_scaled.print();
        // print these: G_bar, R_bar, H_bar;
        // print this: x_cont_init_scaled;
    }
    // TODO pass in gen as a param
    std::vector<i128> sample_knowledge_exponents(i128 from, i128 to_inclusive) {
        i128 N = 6;
        std::vector<i128> res(N);
        for (size_t i = 0; i < N; i++) {
            // TODO update range
            res.at(i) = random_i128(from, to_inclusive);
        }
        // res <- {alpha_0, alpha_1, gamma_0, gamma_1, rho_0, rho_1}
        return res;
    }

    // TODO pass in gen as a param
    std::vector<std::vector<i128>> sample_verification_vectors(i128 m, i128 n, i128 from, i128 to_inclusive) {
        i128 N = 3;
        std::vector<std::vector<i128>> res(N);
        for (size_t i = 0; i < N - 1; i++) {
            res.at(i) = std::vector<i128>(n);
        }
        res.at(N - 1) = std::vector<i128>(m);

        for (auto& x : res) {
            for (size_t i = 0; i < x.size(); i++) {
                x.at(i) = random_i128(from, to_inclusive);
            }
        }
        // res <- {r_0, r_1, s}
        return res;
    }
};

class Encryptor {
private:
i128 v, d, N, q;
std::vector<i128> sk;

public:
    Encryptor(i128 v_, i128 d_, i128 N_, i128 q_, std::vector<i128> sk_)
        : v(v_), d(d_), N(N_), q(q_), sk(sk_) {}
    // Encrypts an RLWE ciphertext of message m
    // m: message polynomial (poly), N: degree, sk: secret key, q: modulus, dth_pows: unused here
    rlwe encrypt_rlwe(const poly& m) {
        poly noise = poly(N);
        auto noise_vec = sample_noise_polynomial(N);
        for (size_t i = 0; i < N; ++i) noise.set(i, noise_vec[i]);

        poly a = poly(N);
        auto a_vec = sample_random_polynomial(N, q);
        for (size_t i = 0; i < N; ++i) a.set(i, a_vec[i]);

        // poly_mult_mod(a, sk, q)
        poly ask(N);
        for (size_t i = 0; i < N; ++i) {
            i128 sum = 0;
            for (size_t j = 0; j < N; ++j) {
                sum += a.get(j) * sk[(i - j + N) % N];
            }
            ask.set(i, mod_(sum, q));
        }

        poly b(N);
        for (size_t i = 0; i < N; ++i) {
            i128 val = m.get(i) + noise.get(i) + ask.get(i);
            b.set(i, mod_(val, q));
        }

        rlwe res(2, N);
        res.set(0, b);
        res.set(1, a);
        return res;
    }

    // Encrypts an RGSW ciphertext of message M (poly)
    rgsw encrypt_rgsw(const poly& M) {
        // Compute v powers: v^0, v^1, ..., v^{d-1}
        std::vector<i128> v_powers(d);
        for (i128 i = 0; i < d; ++i) {
            v_powers[i] = 1;
            for (i128 j = 0; j < i; ++j) {
                v_powers[i] = (v_powers[i] * v) % q;
            }
        }

        // Build G matrix: 2 x 2d, each row is [v_powers, 0...], [0..., v_powers]
        std::vector<std::vector<i128>> G(2, std::vector<i128>(2 * d, 0));
        for (i128 i = 0; i < d; ++i) {
            G[0][i] = v_powers[i];
            G[1][d + i] = v_powers[i];
        }

        // Encryptions of zero
        std::vector<rlwe> encs_of_zero(2 * d);
        poly zero_poly(N);
        for (i128 i = 0; i < 2 * d; ++i) {
            encs_of_zero[i] = encrypt_rlwe(zero_poly);
        }

        // Compute M * G and add to RLWE encryptions
        // NOTE needs modification for
        // FIXME
        for (i128 i = 0; i < N; ++i) {
            for (i128 j = 0; j < 2 * d; ++j) {
                // Add M[i] * G[row][j] to b poly of RLWE
                i128 val = (M.get(i) * G[i < d ? 0 : 1][j]) % q;
                rlwe& ct = encs_of_zero[j];
                poly& b_poly = ct.get_poly(0);
                b_poly.set(i, (b_poly.get(i) + val) % q);
            }
        }

        // Pack RLWE encryptions into RGSW
        rgsw res(2 * d, 2, N);
        for (i128 i = 0; i < 2 * d; ++i) {
            res.set(i, encs_of_zero[i]);
        }
        return res;
    }
    // takes an array2d<i128> and returns an encrypted rgsw_mat
    rgsw_mat encrypt_rgsw_mat(const array2d<i128>& mat) {
        size_t rows = mat.n_rows();
        size_t cols = mat.n_cols();
        rgsw_mat res(rows, cols, 2 * d, 2, N);
        for (size_t i = 0; i < rows; ++i) {
            for (size_t j = 0; j < cols; ++j) {
                poly p(N);
                for (size_t k = 0; k < N; ++k) {
                    p.set(k, mat.get(i, j)); // fill all coeffs with mat(i, j)
                }
                rgsw enc_rgsw = encrypt_rgsw(p);
                res.set(i, j, enc_rgsw);
            }
        }
        return res;
    }

    // Encodes an RGSW plaintext (not encryption, just encoding)
    rgsw encode_rgsw(const poly& M) {
        // Compute v powers: v^0, v^1, ..., v^{d-1}
        std::vector<i128> v_powers(d);
        for (i128 i = 0; i < d; ++i) {
            v_powers[i] = 1;
            for (i128 j = 0; j < i; ++j) {
                v_powers[i] = (v_powers[i] * v) % q;
            }
        }

        // Build G matrix: 2 x 2d, each row is [v_powers, 0...], [0..., v_powers]
        std::vector<std::vector<i128>> G(2, std::vector<i128>(2 * d, 0));
        for (i128 i = 0; i < d; ++i) {
            G[0][i] = v_powers[i];
            G[1][d + i] = v_powers[i];
        }

        // Create rgsw object
        rgsw res(2 * d, 2, N);
        for (i128 j = 0; j < 2 * d; ++j) {
            rlwe ct(2, N);
            for (i128 row = 0; row < 2; ++row) {
                poly p(N);
                for (i128 i = 0; i < N; ++i) {
                    i128 val = (M.get(i) * G[row][j]) % q;
                    p.set(i, val);
                }
                ct.set(row, p);
            }
            res.set(j, ct);
        }
        return res;
    }
};

void test() {
    rlwe_vec u(2, 2, 2);
    for (size_t i = 0; i < u.size(); i++) {
        auto& x = u.get(i);
        for (size_t j = 0; j < x.size(); j++) {
            auto& y = x.get(j);
            for (size_t k = 0; k < y.size(); k++) {
                y.set(k, k);
            }
        }
    }
    std::cout << "u:\n";
    u.print();
    veri_vec_scalar r(2);
    for (size_t i = 0; i < r.size(); i++) {
        r.set(i, i + 1);
    }
    std::cout << "r:\n";
    r.print();
    std::cout << "\n";
    auto ru = r * u;

    std::cout << "ru:\n";
    ru.print();
    i128 expected[] = {0, 3, 6, 9};
    for (auto& x : ru) {
        for (size_t i = 0; i < x.size(); i++) {
            assert(x.get(i) == expected[i]);
        }
    }
    i128 expected1[] = {1, 18, 2, 13};
    rlwe gru = ru.pow();
    std::cout << "gru:\n";
    gru.print();
    for (auto& x : gru) {
        for (size_t i = 0; i < x.size(); i++) {
            assert(x.get(i) == expected1[i]);
        }
    }
    veri_vec_scalar gr = r.pow();
    std::cout << "gr:\n";
    gr.print();
    std::cout << "\n";
    i128 expected2[] = {4, 16};
    for (size_t i = 0; i < r.size(); i++) {
        assert(gr.get(i) == expected2[i]);
    }
    rlwe_vec gu = u.pow();
    std::cout << "gu:\n";
    gu.print();
    i128 expected3[] = {1, 4, 16, 18};
    for (auto& rlwe_el : gu)
        for (auto& poly_el : rlwe_el)
            for (size_t i = 0; i < poly_el.size(); i++) {
                assert(poly_el.get(i) == expected3[i]);
            }

    std::cout << "gru1:\n";
    rlwe gru1 = gr.pow(u);
    gru1.print();

    poly p(2);
    p.set(0, 1);
    p.set(1, 2);
    rlwe_decomp rd(2, 2);
    for (size_t i = 0; i < rd.size(); i++)
        rd.set(i, p); // set all decompositions to the same poly
    std::cout << "rdv:\n";
    rd.print();
    rd.get(0).set(0, 3); // change first coeff of first poly
    std::cout << "rdv after change:\n";
    rd.print();
}

void test_rlwe_decomp() {
    size_t n_coeffs = 2;
    size_t n_polys = 2;
    size_t n_rlwes = 2;
    rlwe_decomp rd(n_polys, n_coeffs);
    for (size_t i = 0; i < rd.size(); i++) {
        poly p(n_coeffs);
        for (size_t j = 0; j < p.size(); j++) {
            p.set(j, j + 1);
        }
        rd.set(i, p);
    }
    std::cout << "rlwe_decomp:\n";
    rd.print();
    rgsw rg(n_rlwes, n_polys, n_coeffs);
    for (size_t i = 0; i < rg.size(); i++) {
        rlwe r(n_polys, n_coeffs);
        for (size_t j = 0; j < r.size(); j++) {
            poly p(n_coeffs);
            for (size_t k = 0; k < p.size(); k++) {
                p.set(k, k + 1);
            }
            r.set(j, p);
        }
        rg.set(i, r);
    }
    std::cout << "rgsw:\n";
    rg.print();

    rlwe rl = rg * rd; // rgsw * rlwe_decomp
    std::cout << "rlwe from rgsw * rlwe_decomp:\n";
    rl.print();

    i128 r = 2;
    rlwe rrl = rl * r; // scalar * rlwe
    std::cout << "rlwe from scalar * rlwe:\n";
    rrl.print();


    rlwe grrl = rrl.pow();
    std::cout << "grrl";
    grrl.print();
    i128 gr = 16; // g^r = 4^2 % 23 = 16
    rlwe grl = rl.pow(gr); // g^r.pow(rlwe)
    std::cout << "grl:\n";
    grl.print();
}

void test_rlwe_decomp_vec() {
    size_t n_coeffs = 2;
    size_t n_polys = 2;
    size_t n_rlwes = 2;
    size_t n_rlwe_decomps = 2;
    rlwe_decomp_vec x(n_rlwe_decomps, n_polys, n_coeffs);
    i128 counter = 0;
    for (size_t v = 0; v < x.size(); v++) {
        rlwe_decomp& rd = x.get(v);
        for (size_t i = 0; i < rd.size(); i++) {
            poly& p = rd.get(i);
            for (size_t j = 0; j < p.size(); j++) {
                p.set(j, counter++ % FIELD_MODULUS); // initialize coefficients to j + 1
            }
            rd.set(i, p);
        }
        x.set(v, rd);
    }
    std::cout << "x:\n";
    x.print();
    size_t rows = 2;
    size_t cols = 2;
    counter = 0;
    rgsw_mat F(rows, cols, n_rlwes, n_polys, n_coeffs);
    for (size_t m = 0; m < F.n_rows(); m++) {
        for (size_t n = 0; n < F.n_cols(); n++) {
            rgsw& rg = F.get(m, n);
            for (size_t i = 0; i < rg.size(); i++) {
                rlwe& r = rg.get(i);
                for (size_t j = 0; j < r.size(); j++) {
                    poly& p = r.get_poly(j);
                    for (size_t k = 0; k < p.size(); k++) {
                        p.set(k, (counter++) % FIELD_MODULUS);
                    }
                }
            }
        }
    }
    std::cout << "F:\n";
    for (size_t i = 0; i < F.n_rows(); i++) {
        for (size_t j = 0; j < F.n_cols(); j++) {
            F.get(i, j).print();
            std::cout << "\n";
        }
    }

    std::cout << "Fx:\n";
    rlwe_vec Fx = F * x; // rgsw_mat * rlwe_decomp_vec
    Fx.print();

    veri_vec_scalar r(n_rlwes);
    for (size_t i = 0; i < r.size(); i++) {
        r.set(i, i + 1);
    }
    std::cout << "r:\n";
    r.print();

    rlwe rFx = r * Fx; // veri_vec_scalar * rlwe_vec
    std::cout << "rFx:\n";
    rFx.print();

    rgsw_vec rF = r * F; // veri_vec_scalar * rgsw_mat
    std::cout << "rF:\n";
    rF.print();


    rgsw_vec grF = rF.pow(); // g^r
    std::cout << "grF:\n";
    grF.print();
    rlwe grFx = grF.pow(x); // g^r.pow(Fx)
    std::cout << "grFx:\n";
    grFx.print();

    rlwe rFx1 = rF * x; // rgsw_vec * rlwe_decomp_vec
    std::cout << "rFx1:\n";
    rFx1.print();

    rlwe grFx1 = rFx1.pow(); // g^rFx1
    std::cout << "grFx1:\n";
    grFx1.print();

    for (size_t i = 0; i < grFx.size(); i++) {
        poly& p1 = grFx.get_poly(i);
        poly& p2 = grFx1.get_poly(i);
        for (size_t j = 0; j < p1.size(); j++) {
            assert(p1.get(j) == p2.get(j)); // check if grFx and grFx1 are equal
        }
    }
}

void test_full() {
    size_t n_rlwes = 7;
    size_t n_polys = 2; // XXX must be 2
    size_t n_coeffs = 4096;
    size_t rows = 5;
    size_t cols = 12;
    rgsw_mat F(rows, cols, n_rlwes, n_polys, n_coeffs);
    rlwe_decomp_vec x(cols, n_rlwes, n_coeffs);
    veri_vec_scalar r(rows);
    init(F, x, r);
    // for (size_t i = 0; i < F.n_rows(); i++) {
    //     for (size_t j = 0; j < F.n_cols(); j++) {
    //         F.get(i, j).print();
    //         std::cout << "\n";
    //     }
    // }
    // x.print();
    // r.print();
    // std::cout << "\n";

    // Compute Fx = F * x
    rlwe_vec Fx = F * x;

    // Compute rF = r * F
    rgsw_vec rF = r * F;

    // Compute grF = g^{rF}
    rgsw_vec grF = rF.pow();

    auto start = std::chrono::high_resolution_clock::now();
    // Compute grFx = (g^{rF})^x = g^{rF * x}
    rlwe grFx = grF.pow(x);
    auto end = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double, std::nano> elapsed = end - start;
    std::cout << "Elapsed time: " << elapsed.count() << " ns\n";
    // init(F, x, r); // reinitialize F, x, r

    // Compute rFx = rF * x
    rlwe rFx = rF * x;
    // rlwe rFx1 = r * Fx;

    // Compute grFx2 = g^{rFx}
    rlwe grFx2 = rFx.pow();

    // init(F, x, r); // reinitialize F, x, r

    // Compute gr = g^r (elementwise)
    veri_vec_scalar gr = r.pow();

    // Compute grFx3 = (g^r)^{Fx} = g^{r * Fx}
    rlwe grFx3 = gr.pow(Fx);

    // std::cout << "grFx:\n";
    // for (const auto& poly_el : grFx) {
    //     for (const auto& coeff : poly_el)
    //         std::cout << print_to_string_i128(coeff) << ", ";
    //     std::cout << "\n";
    // }
    // std::cout << "grFx2:\n";
    // for (const auto& poly_el : grFx2) {
    //     for (const auto& coeff : poly_el)
    //         std::cout << print_to_string_i128(coeff) << ", ";
    //     std::cout << "\n";
    // }
    // std::cout << "grFx3:\n";
    // for (const auto& poly_el : grFx3) {
    //     for (const auto& coeff : poly_el)
    //         std::cout << print_to_string_i128(coeff) << ", ";
    //     std::cout << "\n";
    // }
    // // std::cout << "r:\n";
    // // r.print();
    // // std::cout << "F:\n";
    // // for (size_t i = 0; i < F.n_rows(); i++) {
    // //     for (size_t j = 0; j < F.n_cols(); j++) {
    // //         F.get(i, j).print();
    // //         std::cout << "\n";
    // //     }
    // // }
    // // std::cout << "rF:\n";
    // // rF.print();

    // assert coeffs are equal for grFx, grFx2, grFx3
    for (size_t i = 0; i < grFx.size(); i++) {
        poly& p1 = grFx.get_poly(i);
        poly& p2 = grFx2.get_poly(i);
        poly& p3 = grFx3.get_poly(i);
        for (size_t j = 0; j < p1.size(); j++) {
            assert(p1.get(j) == p2.get(j));
            assert(p1.get(j) == p3.get(j));
        }
    }

}

void run_control_loop() {
    Params pms;
    // params.print();

    using matrix_double = std::vector<std::vector<double>>;
    using matrix_i128 = array2d<i128>;
    using vector_double = std::vector<double>;
    using vector_i128 = std::vector<i128>;
    matrix_double A = pms.A;
    matrix_double B = pms.B;
    matrix_double C = pms.C;
    matrix_i128 F = pms.F;
    matrix_i128 G_bar = pms.G_bar;
    matrix_i128 H_bar = pms.H_bar;
    matrix_i128 R_bar = pms.R_bar;
    vector_double x_plant = pms.x_plant_init;
    poly x_cont = pms.x_cont_init_scaled;

    i128 rr = pms.r;
    i128 ss = pms.s;
    i128 L = pms.L;
    i128 iter_ = pms.iter_;

    i128 from = 1;
    i128 to_inclusive = 100;
    auto knowledge_exps = pms.sample_knowledge_exponents(from, to_inclusive);
    i128 alpha_0 = knowledge_exps.at(0);
    i128 alpha_1 = knowledge_exps.at(1);
    i128 gamma_0 = knowledge_exps.at(2);
    i128 gamma_1 = knowledge_exps.at(3);
    i128 rho_0 = knowledge_exps.at(4);
    i128 rho_1 = knowledge_exps.at(5);
    // std::cout << "Knowledge exponents:\n";
    // for (const auto& exp : knowledge_exps) {
    //     std::cout << print_to_string_i128(exp) << " ";
    // }
    // std::cout << "\n";
    i128 m = 2;
    i128 n = 3;
    auto verification_vectors = pms.sample_verification_vectors(m, n, from, to_inclusive);
    vector_i128 r_0 = verification_vectors.at(0);
    vector_i128 r_1 = verification_vectors.at(1);
    vector_i128 s = verification_vectors.at(2);
    std::cout << "Verification vectors:\n";
    // for (const auto& vec : verification_vectors) {
    //     for (const auto& val : vec) {
    //         std::cout << print_to_string_i128(val) << " ";
    //     }
    //     std::cout << "\n";
    // }
    const i128 N = 4;
    vector_i128 sk = sample_secret_key(N);
    const i128 d = 4;
    double log2q = std::log2(static_cast<double>(FIELD_MODULUS));
    int power = static_cast<int>(std::ceil(log2q / static_cast<double>(d)));
    i128 v = static_cast<i128>(1) << power;

    Encryptor enc(v, d, N, pms.q, sk);
    i128 n_rlwes = 2 * d;
    i128 n_polys = 2;
    i128 n_coeffs = N;
    rgsw_mat F_ctx(F.n_rows(), F.n_cols(), n_rlwes, n_polys, n_coeffs);
    // poly p(n_coeffs);
    // for (size_t i = 0; i < n_coeffs; i++)
    //     p.set(i, 42);
    // rgsw test = enc.encrypt_rgsw(p);
    F_ctx = enc.encrypt_rgsw_mat(F);
}

int main() {
    // test();
    omp_set_nested(1);
    // omp_set_num_threads(1);
    // std::cout << "omp num threads: " << omp_get_max_threads() << std::endl;
    // std::cout << " omp num threads: " << omp_get_num_threads() << std::endl;
    // test_full();


    return 0;
}

// DONE update either i128 to u128, or replace % with mod()
