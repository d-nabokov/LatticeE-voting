from sage.all import *


def _ceil_float(x, pos):
    factor = 10**pos
    return ceil(x * factor) / factor


class PublicParams:
    def __init__(self, Na, Nc, Nv_max):
        self.Na = Na
        self.Nc = Nc
        self.Nv_max = Nv_max
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
        self.kappa = 11
        self.commit_len = self.npoly + 2
        self.baselen = self.lamb + self.kappa + self.commit_len
        self.u = 30

        alpha = 11 * self.k
        self.average_rejection_tries1 = e**(12 / alpha + 1 / (2 * alpha**2))
        theta_1 = _ceil_float(sqrt(101 / (2 * self.baselen * self.d**2 * log(e, 2))) + 5/6, 3)
        self.sigma1 = alpha * sqrt(theta_1 * self.baselen * self.d**2 * self.Na)
        self.beta1 = self.sigma1 * sqrt(2 * self.baselen * self.d)
        self.inf_bound1 = 2**(ceil(log(6 * self.sigma1, 2))) - 1

        self.seedlen = 32

        self.amo_sec_param = 128
        self.amo_n = math.ceil((self.amo_sec_param + 2) / math.log(2 * self.d + 1))
        p = self.number_of_authority_commitments(Nv_max)

        self.average_rejection_tries2 = 3
        theta_2 = sqrt(101 / (2 * self.baselen * self.amo_n * self.d * log(e, 2))) + (5/6) * ((self.u + 1)/self.u)
        theta_2 = _ceil_float(theta_2, 3)
        self.sigma2 = 11 * sqrt(theta_2 * self.baselen * self.amo_n * self.d * p)
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


if __name__ == '__main__':
    PP = PublicParams(2, 5, 10)