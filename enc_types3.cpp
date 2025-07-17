#include <cstdlib>
#include <stdexcept>

template<typename T, typename Derived>
class array1d {
private:
    size_t size_;
    T* arr;
public:
    array1d() : size_(0) {
    }
    array1d(size_t size) : size_(size) {
        arr = new T[size]();
    }
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
    void set(int n, T val) {
        check_index_bounds(n);
        arr[n] = val;
    }
};

class poly : public array1d<int, poly> {
private:
public:
    poly() : array1d<int, poly>() {}
    poly(size_t N) : array1d<int, poly>(N) {
        for (size_t i = 0; i < N; i++) {
            set(i, 0); // initialize coefficients to zero
        }
    }
};

class rlwe : public array1d<poly, rlwe> {
private:
    size_t n_coeffs_;
public:
    rlwe() : array1d<poly, rlwe>() {}
    rlwe(size_t n_polys) : array1d<poly, rlwe>(n_polys) {}
    rlwe(size_t n_polys, size_t n_coeffs) : array1d<poly, rlwe>(n_polys) {
        n_coeffs_ = n_coeffs;
        for (size_t i = 0; i < n_polys; i++) {
            set(i, poly(n_coeffs)); // initialize coefficients to zero
        }
    }
};

class rlwe_vec : public array1d<rlwe, rlwe_vec> {
private:
public:
    rlwe_vec(size_t n_rlwes, size_t n_polys, size_t n_coeffs) : array1d<rlwe, rlwe_vec>(n_rlwes) {
        for (size_t i = 0; i < n_rlwes; i++) {
            set(i, rlwe(n_polys, n_coeffs)); // initialize coefficients to zero
        }
    }
};

int main() {
    rlwe_vec rv(2, 3, 4);
    rv.get(1).get(2).get(3);
    rlwe rl;
    rl.set(0, poly(2));
    return 0;
}