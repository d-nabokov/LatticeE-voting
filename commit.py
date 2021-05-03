from random_polynomials import chi_poly
from linear_alg import matrix_vector, scalar


def commit(PP, B0, b, m, r_seed, nonce):
    r, nonce = chi_poly(PP, PP.baselen, r_seed, nonce)
    t0, t1 = commit_with_r(PP, B0, b, m, r)
    return t0, t1, r, nonce


def commit_with_r(PP, B0, b, m, r):
    t0 = matrix_vector(B0, r, PP.X, PP.d)
    t1 = []
    for i in range(PP.npoly):
        t1.append(scalar(b[i], r, PP.X, PP.d) + m[i])
    return t0, t1