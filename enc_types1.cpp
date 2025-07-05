// TODO how will objects be instantiated setting sizes cascading?
    // perhaps Derived d(size, size_base)
// TODO
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

typedef __int128_t i128;
// global variable for modulus
constexpr i128 FIELD_MODULUS = 23;
constexpr i128 GROUP_MODULUS = 11;
constexpr i128 GENERATOR = 4;

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
                    "(array1d) Value out of range: " + print_to_string_i128(val)
                );
            }
            // else
                // std::cout << "check_value_bounds not implemented for type T\n";
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
                val %= FIELD_MODULUS;
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
            auto val = get(i) * scalar;
            if constexpr (std::is_same_v<T, i128>)
                val %= FIELD_MODULUS;
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
        power %= GROUP_MODULUS;
        base %= FIELD_MODULUS;
        i128 result = 1;
        i128 one = 1;
        i128 two = 2;
        while (power > 0) {
            if ((power % two) == one)
                result = (result * base) % FIELD_MODULUS;
            power >>= 1;
            base = (base * base) % FIELD_MODULUS;
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
            // else
                // std::cout << "check_value_bounds not implemented for type T\n";
        }
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

// TODO constructor to set polys to 1
class poly : public array1d<i128, poly> {
private:
public:
    static constexpr size_t DEFAULT_N = 4;
    poly(size_t N) : array1d<i128, poly>(N) {
        for (size_t i = 0; i < N; i++) {
            // set(i, i + 1); // initialize coefficients to i + 1
            set(i, 0); // initialize coefficients to zero
        }
    }
    poly() : array1d<i128, poly>(DEFAULT_N) {
        for (size_t i = 0; i < DEFAULT_N; i++) {
            // set(i, i + 1); // initialize coefficients to i + 1
            set(i, 0); // initialize coefficients to zero
        }
    }
    auto& get_coeff(size_t n) const {
        return get(n);
    }
};

/* a tuple of polys */
class rlwe : public array1d<poly, rlwe> {
private:
public:
    // NOTE rlwe always has two polys
    static constexpr size_t DEFAULT_N = 2;
    rlwe(size_t N) : array1d<poly, rlwe>(N) {
        assert(N == 2); // rlwe should always have two polys
    }
    rlwe() : array1d<poly, rlwe>(DEFAULT_N) {}
    auto& get_poly(size_t n) const {
        return get(n);
    }
};

class rlwe_vec : public array1d<rlwe, rlwe_vec> {
private:
public:
    static constexpr size_t DEFAULT_N = 3;
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
    static constexpr size_t DEFAULT_N = 4;
    rlwe_decomp_vec(size_t N) : array1d<rlwe_decomp, rlwe_decomp_vec>(N) {}
    rlwe_decomp_vec() : array1d<rlwe_decomp, rlwe_decomp_vec>(DEFAULT_N) {}
    ~rlwe_decomp_vec() {}
};


class rgsw : public array1d<rlwe, rgsw> {
private:
public:
    static constexpr size_t DEFAULT_N = 4;
    rgsw(size_t N) : array1d<rlwe, rgsw>(N) {}
    rgsw() : array1d<rlwe, rgsw>(DEFAULT_N) {}
    ~rgsw() {}

    // NOTE required to expose base class member funcs.
    // Apparently overloads hide them...
    using array1d<rlwe, rgsw>::operator*;
    rlwe operator*(const rlwe_decomp& other) const {
        size_t N = size();
        assert(N == other.size());
        rlwe res(size());
        poly p0(get(0).size());
        poly p1(get(1).size());
        // Initialize all coefficients to zero
        for (size_t j = 0; j < p0.size(); ++j) p0.set(j, 0);
        for (size_t j = 0; j < p1.size(); ++j) p1.set(j, 0);
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
    using array1d<rlwe, rgsw>::pow;

    // Overload for pow with rlwe_decomp argument
    rlwe pow(const rlwe_decomp& other) const {
        size_t N = size();
        assert(N == other.size());
        rlwe res;
        // Initialize p0 and p1 to multiplicative identity (all coefficients 1)
        poly p0(res.get(0).size());
        poly p1(res.get(1).size());
        for (size_t k = 0; k < p0.size(); ++k) {
            p0.set(k, 1);
            p1.set(k, 1);
        }
        for (size_t i = 0; i < N; i++) {
            auto val0 = get(i).get(0).pow(other.get(i));
            auto val1 = get(i).get(1).pow(other.get(i));
            p0 = p0 * val0;
            p1 = p1 * val1;
        }
        res.set(0, p0);
        res.set(1, p1);
        return res;
    }
};

class rgsw_vec : public array1d<rgsw, rgsw_vec> {
private:
public:
    static constexpr size_t DEFAULT_N = 4;
    rgsw_vec(size_t N) : array1d<rgsw, rgsw_vec>(N) {}
    rgsw_vec() : array1d<rgsw, rgsw_vec>(DEFAULT_N) {}
    ~rgsw_vec() {}
    // TODO - [ ] rgsw_vec * rlwe_decomp_vec -> rlwe
    // using array1d<rgsw, rgsw_vec>::operator*;
    using array1d<rgsw, rgsw_vec>::pow;
    rlwe operator*(const rlwe_decomp_vec& other) const {
        size_t n = other.size();
        assert(size() == n);
        rlwe sum; // Initialize sum to default rlwe (all elements zero)
        sum.set(0, poly(other.get(0).size()));
        sum.set(1, poly(other.get(0).size()));
        for (size_t i = 0; i < n; i++) {
            auto val = get(i) * other.get(i);
            sum = sum + val;
        }
        return sum;
    }
    // TODO - [ ] rgsw_vec.pow(rlwe_decomp_vec)->rlwe (needs group mult)
    rlwe pow(const rlwe_decomp_vec& other) const {
        size_t n = other.size();
        assert(size() == n);
        rlwe sum;
        // Initialize sum to multiplicative identity (all polys set to 1)
        for (size_t j = 0; j < sum.size(); ++j) {
            poly ident_poly(sum.get(j).size());
            for (size_t k = 0; k < ident_poly.size(); ++k) {
                ident_poly.set(k, 1);
            }
            sum.set(j, ident_poly);
        }
        for (size_t i = 0; i < n; i++) {
            auto val = get(i).pow(other.get(i));
            sum = sum * val;
        }
        return sum;
    }
};

template<typename T, typename T1, typename U>
U vec_mat_mult(const T& self, const T1& other) {
    size_t rows = self.size();
    size_t cols = self.size(0);
    assert(cols == other.size());
    U res(rows);
    // NOTE get type of array element. Remove reference from type
    using U_element = std::remove_reference_t<decltype(res.get(0))>;
    for (size_t i = 0; i < rows; i++) {
        size_t n_coeffs = self.get(0, 0).get(0).get(0).size();
        size_t n_coeffs_other = other.get(0).get(0).size();
        assert(n_coeffs == n_coeffs_other);
        U_element sum; // Reset sum for each row
        // Initialize sum to additive identity (all polys set to 0)
        for (size_t j = 0; j < n_coeffs; ++j) {
            poly p(n_coeffs);
            sum.set(j, p);
        }
        // for each element in the row
        for (size_t j = 0; j < cols; j++) {
            auto val = self.get(i, j) * other.get(j);
            sum = sum + val;
        }
        res.set(i, sum);
    }
    return res;
}

class rgsw_mat : public array2d<rgsw> {
private:
public:
    constexpr static size_t DEFAULT_M = 3;
    constexpr static size_t DEFAULT_N = 4;
    rgsw_mat(size_t rows, size_t cols) : array2d<rgsw>(rows, cols) {}
    rgsw_mat() : array2d<rgsw>(DEFAULT_M, DEFAULT_N) {}
    ~rgsw_mat() {}
    rlwe_vec operator*(const rlwe_decomp_vec& other) const {
        return vec_mat_mult<rgsw_mat, rlwe_decomp_vec, rlwe_vec>(*this, other);
    }
    rlwe_vec pow(const rlwe_decomp_vec& other) const {
        size_t rows = size();
        size_t cols = size(0);
        assert(cols == other.size());
        rlwe_vec res;
        for (size_t i = 0; i < rows; i++) {
            rlwe sum; // Reset sum for each row
            for (size_t j = 0; j < cols; j++) {
                auto val = get(i, j).pow(other.get(j));
                sum = sum * val;
            }
            res.set(i, sum);
        }
        return res;
    }
};


class veri_vec_scalar : public array1d<i128, veri_vec_scalar> {
private:
public:
    static constexpr size_t DEFAULT_N = 3;
    veri_vec_scalar(size_t N) : array1d<i128, veri_vec_scalar>(N) {}
    veri_vec_scalar() : array1d<i128, veri_vec_scalar>(DEFAULT_N) {}
    ~veri_vec_scalar() {}
    // TODO mult for veri_vec_scalar * rgsw_mat
    // TODO pow func for g^rFx, etc
    rgsw_vec operator*(const rgsw_mat& other) const {
        size_t rows = other.size();
        size_t cols = other.size(0);
        assert(rows == size());
        size_t n_coeffs = other.get(0, 0).get(0).size(); // Get number of coefficients in the first rgsw
        rgsw_vec res(cols);
        for (size_t j = 0; j < cols; j++) {
            rgsw sum(other.get(0, j).size()); // Reset sum for each column
            for (auto& rl : sum) {
                poly p0(n_coeffs);
                poly p1(n_coeffs);
                rl.set(0, p0);
                rl.set(1, p1);
            }
            for (size_t i = 0; i < rows; i++) {
            // for each element in the row
                // scalar mult member function from rgsw_vec
                auto val = other.get(i, j) * get(i);
                sum = sum + val;
            }
            res.set(j, sum);
        }
        return res;
    }
    // TODO - [ ] veri_vec_scalar * rlwe_vec -> rlwe
    rlwe operator*(const rlwe_vec& other) const {
        size_t n = other.size();
        assert(size() == n);
        size_t n_coeffs = other.get(0).get(0).size(); // Get number of coefficients in the first rlwe
        assert(n_coeffs == other.get(0).get(1).size()); // Ensure
        rlwe sum; // Initialize sum to default rlwe (all elements zero)
        // Initialize sum to additive identity (all polys set to 0)
        for (size_t j = 0; j < sum.size(); ++j) {
            poly p(n_coeffs);
            sum.set(j, p);
        }
        for (size_t i = 0; i < n; i++) {
            auto val = other.get(i) * get(i);
            sum = sum + val;
        }
        return sum;
    }
    // TODO - [ ] veri_vec_scalar.pow(rlwe_vec)->rlwe
    using array1d<i128, veri_vec_scalar>::pow;
    rlwe pow(const rlwe_vec& other) const {
        size_t n = other.size();
        assert(size() == n);
        rlwe product;
        // Initialize product to multiplicative identity (all polys set to 1)
        for (size_t j = 0; j < product.size(); ++j) {
            poly ident_poly(product.get(j).size());
            for (size_t k = 0; k < ident_poly.size(); ++k) {
                ident_poly.set(k, 1);
            }
            product.set(j, ident_poly);
        }
        for (size_t i = 0; i < n; i++) {
            auto val = other.get(i).pow(get(i));
            product = product * val;
        }
        return product;
    }

    // rlwe_vec pow(const rlwe_decomp_vec& other) const {
    //     size_t rows = size();
    //     size_t cols = size(0);
    //     assert(cols == other.size());
    //     rlwe_vec res;
    //     rlwe sum;
    //     for (size_t i = 0; i < rows; i++) {
    //         // for each element in the row
    //         for (size_t j = 0; j < cols; j++) {
    //             auto val = get(i, j).pow(other.get(j));
    //             sum = sum * val;
    //         }
    //         res.set(i, sum);
    //     }
    //     return res;
    // }
};

void init(rgsw_mat& F, rlwe_decomp_vec& x, veri_vec_scalar& r) {
    i128 counter = 0;
    for (size_t i = 0; i < F.size(); ++i) {
        for (size_t j = 0; j < F.size(0); ++j) {
            rgsw& rg = F.get(i, j);
            for (size_t k = 0; k < rg.size(); ++k) {
                rlwe& rl = rg.get(k);
                for (size_t l = 0; l < rl.size(); ++l) {
                    poly& p = rl.get_poly(l);
                    for (size_t m = 0; m < p.size(); ++m) {
                        p.set(m, counter++ % FIELD_MODULUS); // just an example initialization
                    }
                }
            }
        }
    }
    // init coeffs of x to i = 0, i++
    counter = 0;
    for (size_t i = 0; i < x.size(); ++i) {
        rlwe_decomp& rd = x.get(i);
        for (size_t j = 0; j < rd.size(); ++j) {
            poly& p = rd.get(j);
            for (size_t m = 0; m < p.size(); ++m) {
                p.set(m, counter++ % FIELD_MODULUS); // just an example initialization
            }
        }
    }
    // init r to i = 0, i++
    counter = 0;
    for (size_t i = 0; i < r.size(); ++i) {
        r.set(i, counter++ % FIELD_MODULUS); // just an example initialization
    }
}

void print_rgsw_mat() {

}

void test() {
    rlwe_vec u(2);
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
    rlwe_decomp rd(2);
    for (size_t i = 0; i < rd.size(); i++)
        rd.set(i, p); // set all decompositions to the same poly
    std::cout << "rdv:\n";
    rd.print();
    rd.get(0).set(0, 3); // change first coeff of first poly
    std::cout << "rdv after change:\n";
    rd.print();
}

void test_rlwe_decomp() {
    rlwe_decomp rd(2);
    for (size_t i = 0; i < rd.size(); i++) {
        poly p(2);
        for (size_t j = 0; j < p.size(); j++) {
            p.set(j, j + 1); // initialize coefficients to j + 1
        }
        rd.set(i, p);
    }
    std::cout << "rlwe_decomp:\n";
    rd.print();
    // rlwe_decomp rd1 = rd.pow(rd);
    // std::cout << "rlwe_decomp after pow:\n";
    // rd1.print();
    rgsw rg(2);
    for (size_t i = 0; i < rg.size(); i++) {
        rlwe r(2);
        for (size_t j = 0; j < r.size(); j++) {
            poly p(2);
            for (size_t k = 0; k < p.size(); k++) {
                p.set(k, k + 1); // initialize coefficients to k + 1
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
}

void test_rlwe_decomp_vec() {
    rlwe_decomp_vec x(2);
    for (size_t v = 0; v < x.size(); v++) {
        rlwe_decomp rd(2);
        for (size_t i = 0; i < rd.size(); i++) {
            poly p(2);
            for (size_t j = 0; j < p.size(); j++) {
                p.set(j, j + 1); // initialize coefficients to j + 1
            }
            rd.set(i, p);
        }
        x.set(v, rd);
    }
    std::cout << "x:\n";
    x.print();
    // rlwe_decomp rd1 = rd.pow(rd);
    // std::cout << "rlwe_decomp after pow:\n";
    // rd1.print();
    rgsw_mat F(2, 2); // 2x2 matrix of rgsw
    for (size_t m = 0; m < F.size(); m++) {
        for (size_t n = 0; n < F.size(0); n++) {
            rgsw rg(2);
            for (size_t i = 0; i < rg.size(); i++) {
                rlwe r(2);
                for (size_t j = 0; j < r.size(); j++) {
                    poly p(2);
                    for (size_t k = 0; k < p.size(); k++) {
                        p.set(k, k + 1); // initialize coefficients to k + 1
                    }
                    r.set(j, p);
                }
                rg.set(i, r);
            }
            F.set(m, n, rg);
        }
    }
    std::cout << "F:\n";
    for (size_t i = 0; i < F.size(); i++) {
        for (size_t j = 0; j < F.size(0); j++) {
            F.get(i, j).print();
            std::cout << "\n";
        }
    }

    std::cout << "Fx:\n";
    rlwe_vec Fx = F * x; // rgsw_mat * rlwe_decomp_vec
    Fx.print();

    veri_vec_scalar r(2);
    for (size_t i = 0; i < r.size(); i++) {
        r.set(i, i + 1); // initialize coefficients to i + 1
    }
    std::cout << "r:\n";
    r.print();

    rlwe rFx = r * Fx; // veri_vec_scalar * rlwe_vec
    std::cout << "rFx:\n";
    rFx.print();

    rgsw_vec rF = r * F; // veri_vec_scalar * rgsw_mat
    std::cout << "rF:\n";
    rF.print();

    rlwe rFx1 = rF * x; // rgsw_vec * rlwe_decomp_vec
    std::cout << "rFx1:\n";
    rFx1.print();
}

void test_full() {
    rgsw_mat F;
    rlwe_decomp_vec x;
    veri_vec_scalar r;
    init(F, x, r);
    for (size_t i = 0; i < F.size(); i++) {
        for (size_t j = 0; j < F.size(0); j++) {
            F.get(i, j).print();
            std::cout << "\n";
        }
    }
    x.print();
    r.print();
    std::cout << "\n";
    // Compute Fx = F * x
    rlwe_vec Fx = F * x;

    // Compute rF = r * F
    rgsw_vec rF = r * F;

    // Compute grF = g^{rF}
    rgsw_vec grF = rF.pow();

    // Compute grFx = (g^{rF})^x = g^{rF * x}
    rlwe grFx = grF.pow(x);

    init(F, x, r); // reinitialize F, x, r

    // Compute rFx = rF * x
    rlwe rFx = rF * x;
    rlwe rFx1 = r * Fx;

    // Compute grFx2 = g^{rFx}
    rlwe grFx2 = rFx.pow();

    init(F, x, r); // reinitialize F, x, r

    // Compute gr = g^r (elementwise)
    veri_vec_scalar gr = r.pow();

    // Compute grFx3 = (g^r)^{Fx} = g^{r * Fx}
    rlwe grFx3 = gr.pow(Fx);

    std::cout << "grFx:\n";
    for (const auto& poly_el : grFx) {
        for (const auto& coeff : poly_el)
            std::cout << print_to_string_i128(coeff) << ", ";
        std::cout << "\n";
    }
    std::cout << "grFx2:\n";
    for (const auto& poly_el : grFx2) {
        for (const auto& coeff : poly_el)
            std::cout << print_to_string_i128(coeff) << ", ";
        std::cout << "\n";
    }
    std::cout << "grFx3:\n";
    for (const auto& poly_el : grFx3) {
        for (const auto& coeff : poly_el)
            std::cout << print_to_string_i128(coeff) << ", ";
        std::cout << "\n";
    }
    std::cout << "r:\n";
    r.print();
    std::cout << "F:\n";
    for (size_t i = 0; i < F.size(); i++) {
        for (size_t j = 0; j < F.size(0); j++) {
            F.get(i, j).print();
            std::cout << "\n";
        }
    }
    std::cout << "rF:\n";
    rF.print();
}

int main() {
    // test();
    test_full();
    // test_rlwe_decomp();
    // test_rlwe_decomp_vec();
    return 0;
}