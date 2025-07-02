<!-- create 3 x 3 table -->
| Stage | Variables | Description |
|-----------|-------------|-------------|
| Setup | \(rx_{k+1} = r_F x_k + r_G y_k \\ su_k = s_H x_k\) | \(r_F, r_G, s_H\): veri_vec_scalar * rgsw_mat -> rgsw_vec
| Setup | \(g^{rx_{k+1}} = g^{r_F x_k} \cdot g^{r_G y_k} \\ g^{su_k} = g^{s_H x_k}\) | \(g^r, g^s, g^{r_F}, g^{s_H}\) - base class pow (veri_vec_scalar.pow(), rgsw_vec.pow())
| Compute | \(x_{k+1} = Fx_k + Gy_k \\ u_k = Hx_k\) | rgsw_mat * rlwe_decomp_vec, rlwe_vec + rlwe_vec
| Compute | \((g^{r})^{x_{k+1}}, (g^{r_F})^{x_k}, (g^{s_H})^{x_k}\) | veri_vec_scalar.pow(rlwe_vec)->rlwe, rgsw_vec.pow(rlwe_decomp_vec)->rlwe (needs to accumulate)
| Verify | \(su_k, r_G y_k \\ g^{s u_k}, g^{r_G y_k}\) | veri_vec_scalar * rlwe_vec, rgsw_vec * rlwe_vec \(\\\)  rlwe_vec.pow()

##### TODO
- [ ] simplify pow1()
- [ ] operator+: confirm this doesn't accumulate
- [ ] should accumulate:
    - [ ] veri_vec_scalar * rgsw_mat -> rgsw_vec

    - [ ] veri_vec_scalar * rlwe_vec -> rlwe
    - [ ] veri_vec_scalar.pow(rlwe_vec)->rlwe

    - [ ] rgsw_vec * rlwe_decomp_vec -> rlwe
    - [ ] rgsw_vec.pow(rlwe_decomp_vec)->rlwe (needs group mult)

    - [ ] rgsw_mat * rlwe_decomp_vec -> rlwe_vec
- [ ] this.pow()
- [ ] rgsw_vec
- [ ] this.pow() should take underying poly elements and do \(g^{element}\)
- [ ] rgsw_vec.pow(rlwe_decomp_vec)
- [ ] rgsw_vec * (rlwe_decomp_vec)
- [ ] group multiplication
- [ ] standard mult requires (mod q)

###### TODO later
- [ ] decomp func
- [ ] PRNG
- [ ] enc/dec
