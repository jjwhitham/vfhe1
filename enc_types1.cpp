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

typedef __int128_t i128;
// constexpr i128 GROUP_MODULUS = 23;
// constexpr i128 FIELD_MODULUS = 11;
// constexpr i128 GENERATOR = 4;
constexpr i128 GROUP_MODULUS = 540431955285196831;
constexpr i128 FIELD_MODULUS = 18014398509506561;
constexpr i128 GENERATOR = 1073741824;

std::string print_to_string_i128(i128 n) {
    if (n == 0) {
        return "0";
    }
    bool neg = false;
    if (n < 0) {
        neg = true;
        n = -n;
    }
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
                    "(array1d) Value out of range: " + print_to_string_i128(val)
                );
            }
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
        check_index_bounds(row, col);
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
        check_index_bounds(row, col);
        check_value_bounds(val);
        arr[row][col] = val;
    }
    size_t get_rows() const {
        return rows_;
    }
    size_t get_cols() const {
        return cols_;
    }
};

class poly : public array1d<i128, poly> {
private:
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
        rlwe res(n_polys, n_coeffs);
        res.set_coeffs_to_one();
        poly& p0 = res.get(0);
        poly& p1 = res.get(1);
        for (size_t i = 0; i < N; i++) {
            auto val0 = get(i).get(0).pow(other.get(i));
            auto val1 = get(i).get(1).pow(other.get(i));
            p0 = p0.group_mult(val0);
            p1 = p1.group_mult(val1);
        }
        res.set(0, p0);
        res.set(1, p1);
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
        rlwe sum(n_polys, n_coeffs);
        sum.set_coeffs_to_one();
        for (size_t i = 0; i < n; i++) {
            auto val = get(i).pow(other.get(i));
            sum = sum.group_mult(val);
        }
        return sum;
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
        size_t rows = get_rows();
        size_t cols = get_cols();
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
        size_t rows = get_rows();
        size_t cols = get_cols();
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
        size_t rows = other.get_rows();
        size_t cols = other.get_cols();
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
        // Initialize product to multiplicative identity (all polys set to 1)
        for (auto& p : product) {
            for (size_t k = 0; k < p.size(); ++k) {
                p.set(k, 1);
            }
        }
        for (size_t i = 0; i < n; i++) {
            auto val = other.get(i).pow(get(i));
            product = product.group_mult(val);
        }
        return product;
    }
};

void init(rgsw_mat& F, rlwe_decomp_vec& x, veri_vec_scalar& r) {
    std::random_device rd;  // Non-deterministic seed
    std::mt19937 gen(rd()); // Mersenne Twister engine
    std::uniform_int_distribution<i128> distrib(0, FIELD_MODULUS - 1); // Range: 0 to 100
    i128 random_number = distrib(gen);
    i128 counter = random_number;
    for (size_t i = 0; i < F.get_rows(); ++i) {
        for (size_t j = 0; j < F.get_cols(); ++j) {
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

void print_rgsw_mat() {

}

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
    for (size_t m = 0; m < F.get_rows(); m++) {
        for (size_t n = 0; n < F.get_cols(); n++) {
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
    for (size_t i = 0; i < F.get_rows(); i++) {
        for (size_t j = 0; j < F.get_cols(); j++) {
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
    size_t cols = 5;
    rgsw_mat F(rows, cols, n_rlwes, n_polys, n_coeffs);
    rlwe_decomp_vec x(cols, n_rlwes, n_coeffs);
    veri_vec_scalar r(rows);
    init(F, x, r);
    // for (size_t i = 0; i < F.get_rows(); i++) {
    //     for (size_t j = 0; j < F.get_cols(); j++) {
    //         F.get(i, j).print();
    //         std::cout << "\n";
    //     }
    // }
    // x.print();
    // r.print();
    // std::cout << "\n";

    auto start = std::chrono::high_resolution_clock::now();
    // Compute Fx = F * x
    rlwe_vec Fx = F * x;

    // Compute rF = r * F
    rgsw_vec rF = r * F;

    // Compute grF = g^{rF}
    rgsw_vec grF = rF.pow();

    // Compute grFx = (g^{rF})^x = g^{rF * x}
    rlwe grFx = grF.pow(x);

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
    // // for (size_t i = 0; i < F.get_rows(); i++) {
    // //     for (size_t j = 0; j < F.get_cols(); j++) {
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

    auto end = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double, std::nano> elapsed = end - start;
    std::cout << "Elapsed time: " << elapsed.count() << " ns\n";
}

int main() {
    // test();
    test_full();
    // test_rlwe_decomp();
    // test_rlwe_decomp_vec();
    return 0;
}