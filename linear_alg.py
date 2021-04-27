from sage.all import *

from params import d, X
from ring import signed_zq


def scalar(a, b):
    res = 0
    for i in range(len(a)):
        res += a[i] * b[i]
    return res.mod(X**d + 1)


def matrix_vector(A, x):
    n = len(A)
    m = len(x)
    assert(len(A[0]) == m)
    res = [0] * n
    for i in range(n):
        res[i] = scalar(A[i], x)
    return res


def matrix_matrix(A, X):
    m = len(A)
    n = len(A[0])
    l = len(X[0])
    assert(len(X) == m)
    res = []
    for i in range(m):
        row = [0] * l
        for j in range(l):
            for k in range(n):
                row[j] = A[i][k] * X[k][j]
        res.append(row)
    return res


def scalar_matrix(A, B):
    m, n = A.dimensions()
    res = 0
    for i in range(m):
        for j in range(n):
            for coef_a, coef_b in zip(A[i][j], B[i][j]):
                res += signed_zq(coef_a) * signed_zq(coef_b)
    return res


def inf_norm(z):
    res = 0
    for coef in z:
        coef = signed_zq(coef)
        if abs(coef) > res:
            res = abs(coef)
    return res


def inf_norm_vect(Z):
    res = 0
    for val in Z:
        norm = inf_norm(val)
        if norm > res:
            res = norm
    return res


def inf_norm_matr(Z):
    res = 0
    for row in Z:
        for val in row:
            norm = inf_norm(val)
            if norm > res:
                res = norm
    return res


def l2_norm(x):
    res = 0
    for coef in x:
        res += signed_zq(coef)**2
    return sqrt(res)


def l2_norm_vect(Z):
    res = 0
    for val in Z:
        for coef in val:
            res += signed_zq(coef)**2
    return sqrt(res)


def l2_norm_matr(Z):
    res = 0
    for row in Z:
        for val in row:
            for coef in val:
                res += signed_zq(coef)**2
    return sqrt(res)