from sage.all import *


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
    return root


class PublicParams:
    def __init__(self, Na, Nc, Nv_max):
        self.Na = Na
        self.Nc = Nc
        self.d = 128
        if Nc <= 32:
            self.l = 32
        else:
            self.l = 128
        self.npolies = ceil(Nc / 128)
        if self.l == 32:
            self.q = 1071646529
            self.q_bits = 30
            self.zeta = 36442389
        elif self.l == 128:
            self.q = 1071650561
            self.q_bits = 30
            self.zeta = 13579840
        else:
            self.q = get_q(l, Integer('111111111000000000000000000000', 2))
            self.q_bits = len(bin(q)[2:])
            print(f'q = {self.q}\nq_bits = {self.q_bits}')
            self.zeta = get_zeta(l, q)
            print(f'zeta = {self.zeta}\n')
        assert((self.q - 1) % (4 * self.l) == 2 * self.l)
        Zq = IntegerModRing(q)
        self.zeta = Zq(zeta)
        self.lamb = 10
        self.kappa = 10
        self.commit_len = npolies + 2
        self.baselen = self.lamb + self.kappa + self.commit_len
        self.u = 10

        # really naive bound
        self.beta1 = Na * self.d
        self.logdelta1 = Integer(ceil(log(self.d * self.baselen * self.beta1, 2)))
        self.delta1 = 1 << self.logdelta1

        self.seedlen = 32

        self.amo_sec_param = 128
        self.amo_n = math.ceil((self.amo_sec_param + 2) / math.log(2 * self.d + 1))
        p = Nv_max
        u_power = self.u
        for i in range(ceil(log(Nv_max, self.u))):
            p += ceil(Nv_max / u_power)
            u_power *= u
        self.p = p

        self.sigma2 = 11 * sqrt(self.baselen * self.amo_n * self.d * self.p)
        self.average_rejection_tries = 3
        self.beta2 = 11 * self.baselen * self.d * sqrt(2 * self.amo_n * self.p)

        self.beta_commit_infty = 1


    def poly_ring(self):
        Zq = IntegerModRing(self.q)
        return PolynomialRing(Zq, name='X').objgen()


    def poly_quot_ring(self):
        P, X = poly_ring(self)
        return P.quotient(X**d + 1, 'x').objgen()


# Na = 2
# Nc = 5
# Nv_max = 10

# d = 128
# l = 128
# assert(Nc <= l)
# if l == 32:
#     q = 1071646529
#     q_bits = 30
#     Zq = IntegerModRing(q)
#     zeta = 36442389
# elif l == 128:
#     q = 1071650561
#     q_bits = 30
#     Zq = IntegerModRing(q)
#     zeta = 13579840
# else:
#     q = get_q(l, Integer('111111111000000000000000000000', 2))
#     q_bits = len(bin(q)[2:])
#     print(f'q = {q}\nq_bits = {q_bits}')
#     Zq = IntegerModRing(q)
#     zeta = get_zeta(l, q)
#     print(f'zeta = {zeta}\n')
# LAMBDA = 10
# KAPPA = 10
# COMMIT_LEN = 3
# BASELEN = LAMBDA + KAPPA + COMMIT_LEN
# u = 10


# # really naive bound
# BETA1 = Na * d
# LOGDELTA1 = Integer(ceil(log(d * BASELEN * BETA1, 2)))
# DELTA1 = 1 << LOGDELTA1

# SEEDLEN = 32

# SEC_PARAM = 128
# N = math.ceil((SEC_PARAM + 2) / math.log(2*d + 1))
# p = Nv_max
# u_power = u
# for i in range(ceil(log(Nv_max, u))):
#     p += ceil(Nv_max / u_power)
#     u_power *= u


# sigma = 11 * sqrt(BASELEN * N * d * p)
# M = 3
# BETA2 = 11 * BASELEN * d * sqrt(2 * N * p)

# BETA_COMMIT_INFTY = 1

# assert((q - 1) % (4*l) == 2*l)

# P, X = PolynomialRing(Zq, name='X').objgen()
# R, x = P.quotient(X**d + 1, 'x').objgen()

if __name__ == '__main__':
    pass