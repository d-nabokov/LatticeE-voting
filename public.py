from params import KAPPA, BASELEN, COMMIT_LEN, d
from random_polynomials import random_poly


def gen_public_b(seed):
    B0, b = gen_public_b_with_extra(seed)
    # gen_public_b_with_extra should call gen_public_b but whatever
    return B0, b[0]


def gen_public_b_with_extra(seed):
    B0 = []
    nonce = 0
    for i in range(KAPPA):
        row = [0] * BASELEN
        row[i] = 1
        for j in range(KAPPA, BASELEN):
            row[j], nonce = random_poly(seed, nonce, d)
        B0.append(row)
        nonce += BASELEN
    b = []
    for i in range(COMMIT_LEN):
        bi = [0] * BASELEN
        bi[i + KAPPA] = 1
        for j in range(COMMIT_LEN + KAPPA, BASELEN):
            bi[j], nonce = random_poly(seed, nonce, d)
        b.append(bi)
        nonce += BASELEN
    return B0, b