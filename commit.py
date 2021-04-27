from params import BASELEN
from random_polynomials import chi_poly
from linear_alg import matrix_vector, scalar


def commit(B0, b1, m, r_seed, nonce):
    r, nonce = chi_poly(BASELEN, r_seed, nonce)
    t0 = matrix_vector(B0, r)
    t1 = scalar(b1, r) + m
    return t0, t1, r, nonce


def commit_with_r(B0, b1, m, r):
    t0 = matrix_vector(B0, r)
    t1 = scalar(b1, r) + m
    return t0, t1