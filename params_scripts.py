from sage.all import *


def get_q(l, starting=2):
    Pr = Primes()
    q = starting
    while not (q - 1) % (4*l) == 2*l:
        q = Pr.next(q)
    return q


def get_zeta(l, q):
    root = None
    Z = IntegerModRing(q)
    for i in range(2, q):
        val = Z(i)
        val_l = val**l
        if val_l != 1 and val_l**2 == 1:
            root = val
            break
    return Integer(root)


def get_any_zeta(l, q):
    Z = IntegerModRing(q)
    u = Z(1)
    while True:
        v = Z(u.nth_root(2 * l))
        v_l = v**l
        if v_l != 1 and v_l**2 == 1:
            return Integer(v)


def get_ring_params(l):
    q = get_q(l, Integer('111111111000000000000000000000', 2))
    q_bits = len(bin(q)[2:])
    zeta = get_zeta(l, q)
    print(f'q = {self.q}\nzeta = {self.zeta}\n')
    return q, q_bits, zeta


def _bit_size(x):
    return ceil(log(x, 2))


# return ballot size in bits
def ballot_size(PP):
    d = PP.d
    res = 3 * d * PP.q_bits + PP.k * PP.baselen * d * _bit_size(12 * PP.sigma1) + 256 + \
        PP.Na * ((PP.kappa + PP.npoly) * d * PP.q_bits + 2 * PP.baselen * d)
    return res


# return authority size in bits
def authority_size(PP):
    d = PP.d
    p = PP.number_of_authority_commitments(PP.Nv_max)
    res = 2 * (PP.baselen * PP.amo_n * d * _bit_size(12 * PP.sigma2) + 256) + \
        2 * PP.baselen * d + (p - PP.Nv_max) * d * PP.q_bits + PP.Nv_max
    return res


def max_prob(PP):
    q = PP.q
    l = PP.l
    k = PP.k
    zeta = Integer(PP.zeta)
    res = 1/q
    s = 0
    for j in range(zeta**k):
        p = 1
        for i in range(l // k):
            p *= abs(1/2 + 1/2 * cos(2 * pi * j * zeta**(k * i) / q))
        s += p
    res += (2 * l // k) / q * s
    return res


def max_prob_test(l):
    q = 4294962689
    if l == 1:
        zeta = (-1) % q
    else:
        zeta = get_any_zeta(l, q)
    print(f'q = {q}, zeta={zeta}')
    res = 1/q
    s = 0
    for j in range(1, (q-1)//(2*l) + 1):
        p = 1
        for k in range(l):
            p *= abs(1/3 + 2/3 * cos(2 * pi * j * zeta**(k) / q))
        s += p
    res += (2 * l) / q * s
    return res


if __name__ == '__main__':
    from params import PublicParams
    # pr = max_prob_test(32)
    # print(pr.n())
    # print(log(pr, 2).n())
    # PP = PublicParams(2, 128, 10)
    # print(log(max_prob(PP), 2).n())
    Na = 4
    for Nc in (32, 128, 1024):
        print(f'Using Na={Na} and Nc={Nc}')
        maxi = 7
        for i in range(1, maxi):
            Nv_max = 10**i
            PP = PublicParams(Na, Nc, Nv_max)
            bs = ballot_size(PP)
            aus = authority_size(PP)
            print(f'Maximum voters={Nv_max:{maxi}}; ballot size = {(bs / (1024 * 8)).n()} KB, authority size = {(aus / (1024 * 8)).n()} KB')
    