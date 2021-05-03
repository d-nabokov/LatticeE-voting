from sage.misc.prandom import random as random_prob
from sage.all import *

from linear_alg import scalar_matrix_zz, scalar_vector_zz, l2_norm_matr, l2_norm_vect
from ring import INTT
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


# returns 1, if Z is independent from B
def rejection_sampling_matrix(Z, B, sigma, M, q):
    scalar = scalar_matrix_zz(Z, B, q)
    norm = l2_norm_matr(B, q)**2
    return _rejection_sampling(scalar, norm, sigma, M)


# returns 1, if z is independent from cr
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


def m_from_vote_arr(PP, vote_arr):
    assert(len(vote_arr) == PP.Nc)
    vlen = PP.npoly * PP.l
    v = vote_arr + [0] * (vlen - PP.Nc)

    return list(INTT(PP, v[i * PP.l : (i + 1) * PP.l]) for i in range(PP.npoly))