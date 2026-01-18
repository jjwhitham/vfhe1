def convolution(a, b, q):
    n = len(a)
    if n != len(b):
        raise ValueError("Polynomials must have the same degree")
    return np.convolve(a, b) % q
def external_product_convolution(C, x, q):
    # TODO update decomp to return np.arrays nested
    assert(C.shape[0] == x.shape[0] * x.shape[1])
    res = np.zeros((2, 2 * C.shape[2] - 1), dtype=np.object_)
    for i in range(res.shape[0]):
        for j in range(x.shape[0]):
            for k in range(x.shape[1]):
                res[i] = (res[i] + convolution(C[j * x.shape[1] + k][i], x[j][k], q)) % q
    return (res % q).astype(np.object_)
def decompose_polynomial(poly, d=7, q=18014398509506561):
    # v - power of 2, s.t. v^{d-1} < q < v^d
    power = int(np.ceil(np.log2(q) / d))
    v = 2**power
    assert(d >= 1)
    # FIXME not sure if we really need the lower bounds check. Fails sometimes
    # assert(v**(d-1) < q and q <= v**d)
    assert(q <= v**d)
    assert(np.all(poly < q))
    assert(np.all(poly >= 0))
    # decompose poly into d polynomials of degree d-1
    polys = []
    N = len(poly)
    for i in range(d):
        polys.append(np.zeros(N, dtype=np.object_))
    for i in range(N):
        for j in range(d):
            polys[j][i] = (v - 1) & (poly[i] // v**j)
    assert(np.all(np.max(polys, axis=1) < v))
    return np.array(polys).astype(np.object_)

def decomp_rlwe(rlwe, d, q):
    """takes an RLWE element and decomposes it"""
    return np.array([decompose_polynomial(rlwe[0], d, q), decompose_polynomial(rlwe[1], d, q)]).astype(np.object_)

def mat_vec_mult_enc(M, v1, d, q):
    """
    Multiplies an encrypted matrix M by an encrypted vector v1 (mod q).
    M is expected to be a list of lists of encrypted ciphertexts (e.g., RGSW ciphertexts),
    and v1 is a list of encrypted RLWE ciphertexts. The function performs decryption-related
    steps such as inverse NTT and decomposition on v1 before computing the external product.
    """
    if len(M[0]) != len(v1):
        raise ValueError("Vector length must match number of columns in matrix")
    # v1_intt = [np.array([np.array(sp.intt(poly, q), dtype=np.object_) for poly in rlwe], dtype=np.object_) for rlwe in v1]
    # v1_intt = [np.array([np.array(intt_1(poly, q, dth_inv_pows), dtype=np.object_) for poly in rlwe], dtype=np.object_) for rlwe in v1]
    # v = [decomp_rlwe(ctx, d, q) for ctx in v1_intt]
    v = [decomp_rlwe(ctx, d, q) for ctx in v1]
    # TODO check this
    # for decomp in v:
    #     for rlwe in decomp:
    #         for i in range(len(rlwe)):
    #             # rlwe[i] = np.array(sp.ntt(rlwe[i], q), dtype=np.object_)
    #             rlwe[i] = np.array(ntt_1(rlwe[i], q, dth_pows), dtype=np.object_)
    return [sum(external_product_convolution(M[i][j], v[j], q) % q for j in range(len(v))) % q for i in range(len(M))]
