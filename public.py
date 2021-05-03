from random_polynomials import random_poly


def gen_public_b(PP, seed):
    B0, b = gen_public_b_with_extra(PP, seed)
    # gen_public_b_with_extra should call gen_public_b but whatever
    return B0, b[:PP.npoly]


def gen_public_b_with_extra(PP, seed):
    B0 = []
    nonce = 0
    for i in range(PP.kappa):
        row = [0] * PP.baselen
        row[i] = 1
        for j in range(PP.kappa, PP.baselen):
            row[j], nonce = random_poly(PP, seed, nonce, PP.d)
        B0.append(row)
        nonce += PP.baselen
    b = []
    for i in range(PP.commit_len):
        bi = [0] * PP.baselen
        bi[i + PP.kappa] = 1
        for j in range(PP.commit_len + PP.kappa, PP.baselen):
            bi[j], nonce = random_poly(PP, seed, nonce, PP.d)
        b.append(bi)
        nonce += PP.baselen
    return B0, b