from sage.all import *

from params import P, X, zeta, l, d, q

def get_q(l, starting=2):
    Pr = Primes()
    q = starting
    while not (q - 1) % (4*l) == 2*l:
        q = Pr.next(q)
    return q


def get_zeta(l, q):
    root = None
    for val in range(2, q):
        if val**(2*l) % q == 1 and not val**(l) % q == 1:
            root = val
            break
    return Zq(root)


def z_star2(x):
    return range(1, x, 2)


def powers_of_zeta(l, k):
    powers = []
    for j in z_star2(2 * l // k):
        for s in range(k):
            p = (j * (2 * l // k + 1)**s) % (2 * l)
            powers.append(p)
    return powers


def NTT(x):
    return list(x.mod(X**(d//l) - zeta**p) for p in powers)


def INTT(r):
    m = list(X**(d//l) - zeta**p for p in powers)
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
    return P(res).mod(X**d + 1)


def signed_zq(v):
    v = Integer(v)
    if v > (q - 1) // 2:
        v = v - q
    return v


powers = powers_of_zeta(l, 1)

if __name__ == '__main__':
    from params import R
    f = R.random_element()
    f = P(f.list())
    assert(INTT(NTT(f)) == f)