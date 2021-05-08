from sage.all import *


class PublicParams:
    def __init__(self, Na, Nc, Nv_max):
        self.Na = Na
        self.Nc = Nc
        self.d = 128
        if Nc <= 32:
            self.l = 32
            self.k = 1
            self.g_zeros = self.d // self.l
        else:
            self.l = 128
            self.k = 4
            self.g_zeros = self.k
        self.npoly = ceil(Nc / self.d)
        if self.l == 32:
            self.q = 1071646529
            self.q_bits = 30
            self.zeta = 36442389
        elif self.l == 128:
            self.q = 1071650561
            self.q_bits = 30
            self.zeta = 13579840
        assert((self.q - 1) % (4 * self.l) == 2 * self.l)
        Zq = IntegerModRing(self.q)
        self.zeta = Zq(self.zeta)
        P, X = PolynomialRing(Zq, name='X').objgen()
        self.P = P
        self.X = X
        R, x = P.quotient(self.X**self.d + 1, 'x').objgen()
        self.R = R

        self.lamb = 10
        self.kappa = 10
        self.commit_len = self.npoly + 2
        self.baselen = self.lamb + self.kappa + self.commit_len
        self.u = 10

        alpha = 11 * self.k
        self.average_rejection_tries1 = e**(12 / alpha + 1 / (2 * alpha**2))
        self.sigma1 = alpha * sqrt(self.baselen * self.d**2 * self.Na)
        self.beta1 = self.sigma1 * sqrt(2 * self.baselen * self.d)
        self.inf_bound1 = 2**(ceil(log(6 * self.sigma1, 2))) - 1

        self.seedlen = 32

        self.amo_sec_param = 128
        self.amo_n = math.ceil((self.amo_sec_param + 2) / math.log(2 * self.d + 1))
        p = self.number_of_authority_commitments(Nv_max)

        self.average_rejection_tries2 = 3
        self.sigma2 = 11 * sqrt(self.baselen * self.amo_n * self.d * p)
        self.beta2 = self.sigma2 * sqrt(2 * self.baselen * self.d)
        self.inf_bound2 = 2**(ceil(log(6 * self.sigma2, 2))) - 1

        self.beta_commit_infty = 1


    def number_of_authority_commitments(self, n):
        p = n
        u_power = self.u
        for i in range(ceil(log(n, self.u))):
            p += ceil(n / u_power)
            u_power *= self.u
        return p


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
    PP = PublicParams(2, 5, 10)