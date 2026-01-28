<!-- create 3 x 3 table -->
| Stage | Variables | Description |
|-----------|-------------|-------------|
| Setup | \(rx_{k+1} = r_F x_k + r_G y_k \\ su_k = s_H x_k\) | \(r_F, r_G, s_H\): veri_vec_scalar * rgsw_mat -> rgsw_vec
| Setup | \(g^{rx_{k+1}} = g^{r_F x_k} \cdot g^{r_G y_k} \\ g^{su_k} = g^{s_H x_k}\) | \(g^r, g^s, g^{r_F}, g^{s_H}\) - base class pow (veri_vec_scalar.pow(), rgsw_vec.pow())
| Compute | \(x_{k+1} = Fx_k + Gy_k \\ u_k = Hx_k\) | rgsw_mat * rlwe_decomp_vec, rlwe_vec + rlwe_vec
| Compute | \((g^{r})^{x_{k+1}}, (g^{r_F})^{x_k}, (g^{s_H})^{x_k}\) | veri_vec_scalar.pow(rlwe_vec)->rlwe, rgsw_vec.pow(rlwe_decomp_vec)->rlwe (needs to accumulate)
| Verify | \(su_k, r_G y_k \\ g^{s u_k}, g^{r_G y_k}\) | veri_vec_scalar * rlwe_vec, rgsw_vec * rlwe_vec \(\\\)  rlwe_vec.pow()

<!-- create 3 x 3 table -->
| Type | Required in | Mbrfuncs |
|-----------|-------------|-------------|
|  | \(rx_{k+1} = r_F x_k + r_G y_k \\ su_k = s_H x_k\) | \(r_F, r_G, s_H\): veri_vec_scalar * rgsw_mat -> rgsw_vec
| Setup | \(g^{rx_{k+1}} = g^{r_F x_k} \cdot g^{r_G y_k} \\ g^{su_k} = g^{s_H x_k}\) | \(g^r, g^s, g^{r_F}, g^{s_H}\) - base class pow (veri_vec_scalar.pow(), rgsw_vec.pow())
| Compute | \(x_{k+1} = Fx_k + Gy_k \\ u_k = Hx_k\) | rgsw_mat * rlwe_decomp_vec, rlwe_vec + rlwe_vec
| Compute | \((g^{r})^{x_{k+1}}, (g^{r_F})^{x_k}, (g^{s_H})^{x_k}\) | veri_vec_scalar.pow(rlwe_vec)->rlwe, rgsw_vec.pow(rlwe_decomp_vec)->rlwe (needs to accumulate)
| Verify | \(su_k, r_G y_k \\ g^{s u_k}, g^{r_G y_k}\) | veri_vec_scalar * rlwe_vec, rgsw_vec * rlwe_vec \(\\\)  rlwe_vec.pow()


##### TODO
Wed 28 Jan
- [ ] add NTL
  - [ ] enc/dec: keep current behaviour, i.e. sample/store a in coeff form



Mon 26 Jan
- [ ] add msm cache test
- [ ] Add invariant: num_pairings

Tue 20th Jan 2026
XXX
- [x] Get asserts working
  - [x] test ops with test_mcl_ops.cpp
- [x] merge commit, rename branch (inc. upstream)


- [ ] Add invariants to counts (calls to funcs, etc)
  - [ ] Add invariant check for get_hash_sec
- XXX [ ] Implement NTTs with NTL


- [x] Implement pairing-based hashing:
  - [x] init G2 and GT
  - [x] replace get_hash_a().pow() with get_hash().pow() (quick hack)
  - [x] pow can be on the hashed types for now
  - [x] ctrl will compute x_d, then call msm() which will compute the msm for each
  poly, setting its G2 val
  - [x] get_hash_sec will then do the pairing
  - [x] use aliases: hashed_t_rgsw_vec -> hashed_rgsw_vec etc to allow switching between the single group and the pairing group implementations

- XXX [ ] CRT can re-use the decomp representation when decomp coeffs are smaller
  than the CRT primes, i.e. for all q_i, v < q_i.

- [x] find BUG:
  - [x] check N: before/after: mod_cyclo, get_hash_a, pow, get_hash_sec
  - [x] check logic for operator*, dot_prod, encode_rgsw

- [ ] Better array1d class
  - [ ] Use inheritance or composition for array1d?
    - [ ] See how composition would work
  - [ ] add 'invariant_size=0'? then: if (invariant_size) assert(size == invariant_size)
  - [ ] remove CRTP: array1d<T, Derived> -> array1d<T>. Extract mbrfncs that use Derived
    - NOTE compiler errors - need to add specialisations

- [ ] Integrate MLC
  - [x] Integrate MLC in the quickest way possible: Do exps in G1
  - [x] Update existing code to flatten hashed RLWEs (use 2m/2n dim for r_i/s)

- [ ] Write performant code
  - [ ] Re-use x_d
  - [x] Fix redundant NTTs: convert x, y, u_re to eval form once
    - [x] Use more distinctive names: mod_cyclo(), convert_to_eval/eval_form, coeff_form
    - [x] Remove some excessive copying
    - [ ] Remove more copying
  - [ ] Compare MCL against BLST and others
  - [ ] Allocate all variables upfront and mutate them across time
    - [ ] Q: How for NTL & MCL?
  - [ ] Investigate increasing stack size and perform everything on stack
  - [ ] Experiment: Fastest scalars - MCL or NTL?
  - [ ] Want bigz type that mods itself once it reaches a threshold (NTL?)

- [ ] Start new codebase (called vHEctrl? vHE? VEctlr?)
  - [ ] If NTL > MCL: Use NTL ZZ for ints ZZX for polys
  - [ ] Make design decision: vec_ZZX for rlwe and std::vector<vec_ZZX>, vecs of vecs of vec_ZZX etc. to make compound types, or just have vec_ZZX and std::vector<vec_ZZX> as the core 1d/2d arrays of polys, with everything else just aliasing these...

The classes Vec<ZZ_p> (a.k.a., vec_ZZ_p), Mat<ZZ_p> (a.k.a., mat_ZZ_p), and ZZ_pX represent vectors, matrices, and polynomials mod p, and work much the same way as the corresponding classes for ZZ.

###### TODO future
- [ ] Build with CMake
- [ ] Add simple testing (Ctest? Might be from Google)
- [ ] address santisation, valgrind and static analysis
- [ ] Investigate hardware acceleration options
