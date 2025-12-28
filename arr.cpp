#include <cstddef>
#include <iostream>
#include <cassert>

template<typename T>
class arr {
    size_t n;
    T* x;
public:
    arr() : n{0}, x{nullptr} {
        std::cout << "arr() called\n";
    }
    arr(size_t n_) : n{n_}, x{new T[n_]} {
        std::cout << "arr(size_t n) called\n";
    }
    // copy ctor
    arr(const arr& other) : n{other.n}, x{new T[other.n]} {
        std::cout << "copy ctor called\n";
        for (size_t i = 0; i < n; i++)
            x[i] = other.x[i];
    }
    // move ctor
    arr(arr&& other) noexcept : n{other.n}, x{other.x} {
        std::cout << "move ctor called\n";
        other.n = 0;
        other.x = nullptr;
    }
    // dtor
    ~arr() {
        std::cout << "~arr() called\n";
        delete[] x;
    }
    // copy assignment
    arr& operator=(const arr& other) {
        std::cout << "copy assignment called\n";
        if (this != &other) { // TODO see cppcon RO5 when this can be dropped
            if (x == nullptr) { // constructed with arr()
                n = other.n;
                x = new T[n];
            } else // constructed with arr(size_t n)
                assert(n == other.n);
            // delete[] x;
            for (size_t i = 0; i < n; i++) {
                if constexpr (std::is_same_v<T, int>)
                    check(x[i]);
                x[i] = other.x[i];
            }
        }
        return *this;
    }
    // move assignment
    arr& operator=(arr&& other) noexcept { // TODO understand noexcept
        std::cout << "move assignment called\n";
        if (this != &other) {
            if (x == nullptr) // constructed with arr()
                n = other.n;
            else { // constructed with arr(size_t n)
                assert(n == other.n);
                delete[] x;
            }
            x = other.x;
            other.n = 0;
            other.x = nullptr;
        }
        return *this;
    }
    arr& operator+=(const arr& other) {
        std::cout << "operator+= called\n";
        for (size_t i = 0; i < n; i++) {
            x[i] += other.x[i];
            if constexpr (std::is_same_v<T, int>)
                x[i] %= 7;
        }
        return *this;
    }
    arr& operator+=(arr&& other) noexcept {
        std::cout << "(move) operator+= called\n";
        for (size_t i = 0; i < n; i++) {
            x[i] += other.x[i];
            if constexpr (std::is_same_v<T, int>)
                x[i] %= 7;
        }
        return *this;
    }
    // // get type of *this as Derived
    // arr operator+(const arr& other) {
    //     std::cout << "operator+ called\n";
    //     assert(n == other.n);
    //     // get type of *this as Derived
    //     // using Derived = std::remove_reference_t<decltype(*this)>;
    //     arr ret{n};
    //     for (size_t i = 0; i < n; i++) {
    //         T val = static_cast<T>(x[i] + other.x[i]);
    //         if constexpr (std::is_same_v<T, int>)
    //             val %= 7;
    //         // ret.get(i) = val;
    //         ret.set(i, val);
    //     }
    //     return ret;
    // }
    T& get(size_t i) const {
        return x[i];
    }
    T& operator[](size_t i) const {
        return x[i];
    }
    void set(size_t i, const T& val) {
        x[i] = val;
    }
    void set(size_t i, T&& val) noexcept {
        x[i] = std::move(val);
    }
    size_t size() const {
        return n;
    }
    void check(const T& val) const {
        (void)val;
    }

    // Member templated operator+ that returns the derived type.
    template<typename Derived>
    Derived operator+(const Derived& b) const {
        std::cout << "(member) Derived operator+ called\n";
        assert(size() == b.size());
        Derived ret(size());
        for (size_t i = 0; i < size(); i++) {
            ret[i] = get(i);
            ret[i] += b.get(i);
        }
        return ret;
    }

};

// // Non-member templated operator+ that returns the derived type.
// template<typename Derived>
// Derived operator+(const Derived& a, const Derived& b) {
//     std::cout << "Derived operator+ called\n";
//     assert(a.size() == b.size());
//     Derived ret(a.size());
//     for (size_t i = 0; i < a.size(); i++) {
//         // auto tmp = a.get(i);
//         // tmp += b.get(i);
//         // ret.set(i, tmp);

//         ret[i] = a.get(i);
//         ret[i] += b.get(i);

//         // ret.set(i, std::move(a.get(i) + b.get(i))); // calls copy assgn
//         // ret.set(i, a.get(i) + b.get(i)); // calls copy assgn
//         // ret.get(i) = a.get(i) + b.get(i); // calls move assign

//         // ret[i] = a[i] + b[i]; // calls move assign
//     }
//     return ret;
// }

class poly : public arr<int> {
public:
    void speak() {
        std::cout << "hello, world!\n";
    }
};
class rlwe : public arr<poly> { };
// using rlwe = arr<poly>;

using rgsw = arr<rlwe>;

int main() {
    // TODO check that a mbrfunc that returns arr& instead of Derived& can still
    // call a Derived mbrfunc
    poly x{2};
    x.set(0, 2);
    x.set(1, 3);
    poly y{2};
    y.set(0, 3);
    y[1] = 4;


    rlwe z{2};
    z.set(0, x);
    z.set(1, y);
    // rlwe zz = std::move(z); // doesn't call move assignment
    // rlwe zz{2} = z + z;
    // rlwe zz = static_cast<rlwe>(z + z);

    rlwe zz = z + z; // calls copy/move assign
    zz[0][1] = 6;

    // rlwe zz;
    // zz = z + z; // calls move assign

    // zz.get(0) = z.get(0);
    // rlwe zz{std::move(z)};
    // rlwe zz{z};
    // zz += z;
    // std::cout << "z.get(1).get(1): " << zz.get(1).get(1) << "\n";

    // x += y;
    // assert(x.get(0) == 5);
    // assert(x.get(1) == 0);
    // x.speak();


    // DONE see if memory is leaked from nested arr (recursive dtor required?)
    //   DONE create some inheritence
    // DONE what are the implications of returning arr& v arr<T, Derived>& v Derived&?

    return 0;
}