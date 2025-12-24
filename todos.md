<!-- create 3 x 3 table -->
| Stage | Variables | Description |
|-----------|-------------|-------------|
| Setup | \(rx_{k+1} = r_F x_k + r_G y_k \\ su_k = s_H x_k\) | \(r_F, r_G, s_H\): veri_vec_scalar * rgsw_mat -> rgsw_vec
| Setup | \(g^{rx_{k+1}} = g^{r_F x_k} \cdot g^{r_G y_k} \\ g^{su_k} = g^{s_H x_k}\) | \(g^r, g^s, g^{r_F}, g^{s_H}\) - base class pow (veri_vec_scalar.pow(), rgsw_vec.pow())
| Compute | \(x_{k+1} = Fx_k + Gy_k \\ u_k = Hx_k\) | rgsw_mat * rlwe_decomp_vec, rlwe_vec + rlwe_vec
| Compute | \((g^{r})^{x_{k+1}}, (g^{r_F})^{x_k}, (g^{s_H})^{x_k}\) | veri_vec_scalar.pow(rlwe_vec)->rlwe, rgsw_vec.pow(rlwe_decomp_vec)->rlwe (needs to accumulate)
| Verify | \(su_k, r_G y_k \\ g^{s u_k}, g^{r_G y_k}\) | veri_vec_scalar * rlwe_vec, rgsw_vec * rlwe_vec \(\\\)  rlwe_vec.pow()

##### TODO
- [ ] Integrate MLC
  - [ ] Get MLC in the quickest way possible: Do exps in G1
  - [ ] Update existing code to flatten hashed RLWEs (use 2m/2n dim for r_i/s)
    - [ ]

- [ ] Write performant code
  - [ ] Fix redundant NTTs: convert x, y, u_re to eval form once
    - [ ] Use more distinctive names: mod_cyclo(), convert_to_eval/eval_form, coeff_form
    - [ ] Remove excessive copying
  - [ ] Compare MCL against BLST and others
  - [ ] Allocate all variables upfront and mutate them across time
    - [ ] Q: How for NTL & MCL?
  - [ ] Investigate increasing stack size and performing everything on stack
  - [ ] Experiment: Fastest scalars - MCL or NTL?

- [ ] Start new codebase (called vHEctrl? vHE? VEctlr?)
  - [ ] If NTL > MCL: Use NTL ZZ for ints ZZX for polys
  - [ ] Make design decision: vec_ZZX for rlwe and std::vector<vec_ZZX>, vecs of vecs of vec_ZZX etc. to make compound types, or just have vec_ZZX and std::vector<vec_ZZX> as the core 1d/2d arrays of polys, with everything else just aliasing these...

The classes Vec<ZZ_p> (a.k.a., vec_ZZ_p), Mat<ZZ_p> (a.k.a., mat_ZZ_p), and ZZ_pX represent vectors, matrices, and polynomials mod p, and work much the same way as the corresponding classes for ZZ.

###### TODO future
- [ ] Build with CMake
- [ ] Add simple testing (Ctest? Might be from Google)
- [ ] Investigate hardware acceleration options
