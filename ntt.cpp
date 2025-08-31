#include "ntt.hpp"

int main() {
    constexpr u128 INV_2ROU = pow_constexpr(TWO_ROU, FIELD_MOD - 2, FIELD_MOD);
    constexpr u128 INV_2N = pow_constexpr(2 * POLY_SIZE, FIELD_MOD - 2, FIELD_MOD);
    // check n, q
    constexpr bool q_and_N_legal = are_q_and_N_legal(FIELD_MOD, POLY_SIZE);

    // check primitive nth ROU, w
    constexpr bool w_is_legal = is_w_legal(NTH_ROU);
    // check primitive 2nth ROU, psi
    constexpr bool psi_is_legal = is_psi_legal(TWO_ROU);
    // compute primitive nth ROU powers, w^i[]
    constexpr arr_u128 w_pows = get_rou_pows(NTH_ROU);
    constexpr arr_u128 psi_pows = get_rou_pows(TWO_ROU);
    constexpr arr_u128 psi_inv_pows = get_rou_pows(INV_2ROU);
    // check rou_pows
    constexpr bool w_pows_is_legal = is_w_pows_legal(w_pows);
    constexpr bool psi_pows_is_legal = is_psi_pows_legal(psi_pows);
    constexpr bool all_legal = q_and_N_legal && w_is_legal && psi_is_legal \
        && w_pows_is_legal && psi_pows_is_legal;
    static_assert(all_legal);

    arr2n_u128 x{};
    for (size_t i = 0; i < 2 * POLY_SIZE; i++)
        assert(x[i] == 0);
    for (size_t i = 0; i < POLY_SIZE; i++) {
        x[i] = i + 1;
    }
    span_u128 x_span(x);
    ntt_recursive(x_span, TWO_ROU);
    // for (size_t i = 0; i < 2 * POLY_SIZE; i++) {
    //     x[i] = x[i] * x[i] % FIELD_MOD;
    // }
    print_arr(x);
    intt_recursive(x_span, INV_2ROU, INV_2N);
    print_arr(x);

    ntt_iter(x, psi_pows);
    print_arr(x);
    intt_iter(x, psi_inv_pows, INV_2N);
    print_arr(x);
    for (size_t i = 0; i < POLY_SIZE; i++)
        assert(x[i] == i + 1);
    for (size_t i = POLY_SIZE; i < 2 * POLY_SIZE; i++)
        assert(x[i] == 0);
    return 0;
}