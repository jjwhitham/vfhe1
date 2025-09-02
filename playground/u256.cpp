/*
This header provides a custom 256 bit unsigned integer type with features:
  - type is a std::array of 2 * __uint128_t allowing static, auto or heap storage
  - the two limbs are in little-endian order
  - u256_mod inherets from u256 and provides automatic mod reductions via the use
    of overloaded operator member functions
*/
#include <array>
#include <cstdlib>

constexpr size_t N_LIMBS_U256 = 2;
using u128 = __uint128_t;

class u256 {
private:
    std::array<u128, N_LIMBS_U256> _val;
public:
    u256() {
        _val[0] = 0;
        _val[1] = 0;
    }

    u256(u128 hi, u128 low) {
        _val[0] = low;
        _val[1] = hi;
    }

    u256 operator+(const u256& other) const {
        // check for overflow of low limb
        u128 low = _val[0] + other[0];
        // check for overflow
        int carry = (low < _val[0]) ? 1 : 0; // && low < other[0] is implied
        u128 hi = _val[1] + other[1] + carry;
        return u256(hi, low);
    }
    //
    u128 operator[](size_t i) const {
        return _val.at(i);
    }

};

int main() {

}