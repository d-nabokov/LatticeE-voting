from random_polynomials import chi_poly
from linear_alg import matrix_vector, scalar


def commit(PP, B0, b1, m, r_seed, nonce):
    r, nonce = chi_poly(PP.baselen, r_seed, nonce)
    t0 = matrix_vector(B0, r, PP.X, PP.d)
    t1 = scalar(b1, r, PP.X, PP.d) + m
    return t0, t1, r, nonce


def commit_with_r(PP, B0, b1, m, r):
    t0 = matrix_vector(B0, r, PP.X, PP.d)
    t1 = scalar(b1, r, PP.X, PP.d) + m
    return t0, t1