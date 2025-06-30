#include "vfhe.h"

class poly : public array1d {
private:
public:
    poly(size_t N) : array1d(N) {
    }
    ~poly() = default;

    poly operator*(const poly& other) const {
        size_t N = size();
        poly res(N);
        for (size_t i = 0; i < N; i++) {
            res[i] = (*this)[i] * other[i];
        }
        return res;
    }
};

/* a tuple of polys */
class rlwe {
    private:
    public:
    rlwe() {}
    ~rlwe() = default;
};

class rlwe_vec {
    private:
    public:
    rlwe_vec() {}
    ~rlwe_vec() {}
};

class rlwe_decomp {
    private:
    public:
    rlwe_decomp() {}
    ~rlwe_decomp() {}
};

class rlwe_decomp_vec {
    private:
    public:
    rlwe_decomp_vec() {}
    ~rlwe_decomp_vec() {}
};

class rgsw {
    private:
    public:
    rgsw() {}
    ~rgsw() {}
};

class rgsw_mat {
    private:
    public:
    rgsw_mat() {}
    ~rgsw_mat() {}
};

int main() {

    i128 arr1[] = { 1, 2, 3, 4 };
    i128 arr2[] = { 1, 2, 3, 4 };
    poly p1(4);
    poly p2(4);
    for (int i = 0; i < 4; i++) {
        p1[i] = arr1[i];
        p2[i] = arr2[i];
    }
    poly p3 = p1 * p2;
    for (int i = 0; i < 4; i++)
        cout << print_to_string_i128(p3[i]) + ", ";
    cout << "\n";
    return 0;
}