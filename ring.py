from sage.all import *


def z_star2(x):
    return range(1, x, 2)


def powers_of_zeta(zeta, l, k):
    powers = []
    for j in z_star2(2 * l // k):
        for s in range(k):
            p = (j * (2 * l // k + 1)**s) % (2 * l)
            powers.append(p)
    return powers


def NTT(PP, x):
    global powers
    if powers is None:
        powers = powers_of_zeta(PP.zeta, PP.l, PP.k)
    return list(x.mod(PP.X**(PP.d//PP.l) - PP.zeta**p) for p in powers)


def INTT(PP, r):
    global powers
    if powers is None:
        powers = powers_of_zeta(PP.zeta, PP.l, PP.k)
    m = list(PP.X**(PP.d//PP.l) - PP.zeta**p for p in powers)
    P, X = PP.P, PP.X
    M = 1
    for poly in m:
        M *= poly
    res = 0
    for i in range(len(r)):
        Mi = M / m[i]
        R, x = P.quotient(m[i], name='x').objgen()

        xi = r[i] * Mi * P((1/R(Mi)).list())
        xi = P(xi)
        # xi = P(r[i].list()) * Mi * P((1/R(Mi)).list())
        res += xi
    return P(res).mod(X**PP.d + 1)


def signed_zq(v, q):
    v = Integer(v)
    if v > (q - 1) // 2:
        v = v - q
    return v


powers = None


if __name__ == '__main__':
    from params import PublicParams
    PP = PublicParams(2, 5, 10)
    R, x = PP.poly_quot_ring()
    f = R.random_element()
    f = PP.P(f.list())
    assert(INTT(PP, NTT(PP, f)) == f)

    # TODO: fix fully-splitting NTT
    PP = PublicParams(2, 127, 10)
    R, x = PP.poly_quot_ring()
    f = R.random_element()
    f = PP.P(f.list())
    assert(INTT(PP, NTT(PP, f)) == f)