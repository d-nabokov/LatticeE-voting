from sage.misc.prandom import random as random_prob
from sage.all import *

from linear_alg import scalar_matrix, l2_norm_matr
from os import urandom


def poly_to_bytes(p):
    return b''.join(int(x).to_bytes(4, byteorder='big') for x in p.list())


def randombytes(l):
    return urandom(l)


def rejection_sampling(Z, B, sigma, M, q):
    u = random_prob()
    scalar = scalar_matrix(Z, B, q)
    norm = l2_norm_matr(B, q)**2
    border = 1 / M * e**((-2*scalar + norm) / (2 * sigma**2))
    if u > border:
        return 0
    else:
        return 1


def max_with_index(l):
    m = 0
    for i, val in enumerate(l):
        if val > m:
            m = val
            index = i
    return m, index