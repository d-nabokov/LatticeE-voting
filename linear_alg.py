from sage.all import *

from ring import signed_zq


def scalar(a, b, X, d):
    res = 0
    for i in range(len(a)):
        res += a[i] * b[i]
    return res.mod(X**d + 1)


def matrix_vector(A, x, X, d):
    n = len(A)
    m = len(x)
    assert(len(A[0]) == m)
    res = [0] * n
    for i in range(n):
        res[i] = scalar(A[i], x, X, d)
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


def scalar_vector_zz(a, b, q):
    n = len(a)
    res = 0
    for i in range(n):
        for coef_a, coef_b in zip(a[i], b[i]):
            res += signed_zq(coef_a, q) * signed_zq(coef_b, q)
    return res


def scalar_matrix_zz(A, B, q):
    m, n = A.dimensions()
    res = 0
    for i in range(m):
        for j in range(n):
            for coef_a, coef_b in zip(A[i][j], B[i][j]):
                res += signed_zq(coef_a, q) * signed_zq(coef_b, q)
    return res


def inf_norm(z, q):
    res = 0
    for coef in z:
        coef = signed_zq(coef, q)
        if abs(coef) > res:
            res = abs(coef)
    return res


def inf_norm_vect(Z, q):
    res = 0
    for val in Z:
        norm = inf_norm(val, q)
        if norm > res:
            res = norm
    return res


def inf_norm_matr(Z, q):
    res = 0
    for row in Z:
        for val in row:
            norm = inf_norm(val, q)
            if norm > res:
                res = norm
    return res


def l2_norm(x, get_q):
    res = 0
    for coef in x:
        res += signed_zq(coef, q)**2
    return sqrt(res)


def l2_norm_vect(Z, q):
    res = 0
    for val in Z:
        for coef in val:
            res += signed_zq(coef, q)**2
    return sqrt(res)


def l2_norm_matr(Z, q):
    res = 0
    for row in Z:
        for val in row:
            for coef in val:
                res += signed_zq(coef, q)**2
    return sqrt(res)