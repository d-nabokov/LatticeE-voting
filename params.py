from sage.all import *

Na = 2
Nc = 5
Nv_max = 10

d = 128
l = 32
assert(Nc <= l)
if l == 32:
    q = 1071646529
    q_bits = 30
    zeta = 36442389
else:
    q = get_q(l, Integer('111111111000000000000000000000', 2))
    q_bits = len(bin(q)[2:])
    zeta = get_zeta(l, q)
    print(f'q = {q}\nzeta = {zeta}\n')
LAMBDA = 10
KAPPA = 10
COMMIT_LEN = 3
BASELEN = LAMBDA + KAPPA + COMMIT_LEN
u = 10


# really naive bound
BETA1 = Na * d
LOGDELTA1 = Integer(ceil(log(d * BASELEN * BETA1, 2)))
DELTA1 = 1 << LOGDELTA1

SEEDLEN = 32

SEC_PARAM = 128
N = math.ceil((SEC_PARAM + 2) / math.log(2*d + 1))
p = Nv_max
u_power = u
for i in range(ceil(log(Nv_max, u))):
    p += ceil(Nv_max / u_power)
    u_power *= u


sigma = 11 * sqrt(BASELEN * N * d * p)
M = 3
BETA2 = 11 * BASELEN * d * sqrt(2 * N * p)

BETA_COMMIT_INFTY = 1

assert((q - 1) % (4*l) == 2*l)

Zq = IntegerModRing(q)
P, X = PolynomialRing(Zq, name='X').objgen()
R, x = P.quotient(X**d + 1, 'x').objgen()

if __name__ == '__main__':
    print('printing from params')