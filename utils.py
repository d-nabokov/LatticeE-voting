from sage.misc.prandom import random as random_prob
from sage.all import *

from linear_alg import scalar_matrix_zz, scalar_vector_zz, l2_norm_matr, l2_norm_vect
from os import urandom


def poly_to_bytes(p):
    return b''.join(int(x).to_bytes(4, byteorder='big') for x in p.list())


def randombytes(l):
    return urandom(l)


def _rejection_sampling(scalar, norm, sigma, M):
    u = random_prob()
    border = 1 / M * e**((-2*scalar + norm) / (2 * sigma**2))
    if u > border:
        return 0
    else:
        return 1


def rejection_sampling_matrix(Z, B, sigma, M, q):
    scalar = scalar_matrix_zz(Z, B, q)
    norm = l2_norm_matr(B, q)**2
    return _rejection_sampling(scalar, norm, sigma, M)


def rejection_sampling_vector(z, cr, sigma, M, q):
    scalar = scalar_vector_zz(z, cr, q)
    norm = l2_norm_vect(cr, q)**2
    return _rejection_sampling(scalar, norm, sigma, M)


def max_with_index(l):
    m = 0
    for i, val in enumerate(l):
        if val > m:
            m = val
            index = i
    return m, index