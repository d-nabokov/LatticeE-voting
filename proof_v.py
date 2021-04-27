from Crypto.Hash import SHAKE128

from params import KAPPA, SEEDLEN, DELTA1, BETA1, BASELEN, d, l, X, Nc
from ring import INTT
from public import gen_public_b_with_extra
from utils import poly_to_bytes, randombytes
from linear_alg import inf_norm, matrix_vector, scalar
from random_polynomials import challenge, random_poly, random_poly_with_zeros, uniform_poly
from commit import commit


def get_gamma(t0, t1, t2, w):
    shake = SHAKE128.new()
    for i in range(KAPPA):
        shake.update(poly_to_bytes(t0[i]))
    shake.update(poly_to_bytes(t1))
    shake.update(poly_to_bytes(t2))
    for i in range(KAPPA):
        shake.update(poly_to_bytes(w[i]))
    gamma_hash = shake.read(int(SEEDLEN))
    gamma, _ = random_poly(gamma_hash, 0, d//l)
    return gamma, gamma_hash


def get_challenge_hash(gamma_hash, t3, vpp, h, vulp):
    shake = SHAKE128.new()
    shake.update(gamma_hash)
    shake.update(poly_to_bytes(t3))
    shake.update(poly_to_bytes(vpp))
    shake.update(poly_to_bytes(h))
    shake.update(poly_to_bytes(vulp))
    c_hash = shake.read(int(SEEDLEN))
    return c_hash


def get_challenge(c_hash):
    return challenge(c_hash)


def check_z_len(z):
    border = DELTA1 - BETA1
    for i in range(BASELEN):
        if inf_norm(z[i]) >= border:
            return 1
    return 0


def proof_v(t0, t1, r, m, Nc, public_seed):
    try_index = 0
    m_prime = INTT([1] * Nc + [0] * (l - Nc))
    B0, b = gen_public_b_with_extra(public_seed)
    seed = randombytes(SEEDLEN)
    nonce = 0
    g, nonce = random_poly_with_zeros(seed, nonce, d, d//l)
    t2 = scalar(b[1], r) + g
    while True:
        print('iteration', try_index)
        try_index += 1
    
        y, nonce = uniform_poly(BASELEN, seed, nonce)
        w = matrix_vector(B0, y)

        gamma, gamma_hash = get_gamma(t0, t1, t2, w)
        t3 = (scalar(b[2], r) - (2 * m - m_prime) * scalar(b[0], y)).mod(X**d + 1)
        vpp = (scalar(b[2], y) + scalar(b[0], y)**2).mod(X**d + 1)
        intt_factor = gamma * l
        h = (g + intt_factor * m - gamma).mod(X**d + 1)
        vulp = scalar(list((intt_factor * b[0][i] + b[1][i]) for i in range(BASELEN)), y)

        c_hash = get_challenge_hash(gamma_hash, t3, vpp, h, vulp)
        c = get_challenge(c_hash)

        z = [0] * BASELEN
        for i in range(BASELEN):
            z[i] = (y[i] + c * r[i]).mod(X**d + 1)
        if not check_z_len(z):
            break
    return (h, c_hash, z), (t2, t3)


def verify_v(proof, commitment, additional_com, public_seed):
    m_prime = INTT([1] * Nc + [0] * (l - Nc))
    B0, b = gen_public_b_with_extra(public_seed)
    h, c_hash, z = proof
    c = get_challenge(c_hash)
    t0, t1 = commitment
    t2, t3 = additional_com
    if check_z_len(z):
        print('check_z_len')
        return 1
    B0z = matrix_vector(B0, z)
    w = [0] * KAPPA
    for i in range(KAPPA):
        w[i] = (B0z[i] - c * t0[i]).mod(X**d + 1)
    f1 = (scalar(b[0], z) - c * t1).mod(X**d + 1)
    f2 = (scalar(b[0], z) - c * (t1 - m_prime)).mod(X**d + 1)
    f3 = (scalar(b[2], z) - c * t3).mod(X**d + 1)
    vpp = (f1 * f2 + f3).mod(X**d + 1)
    hlist = h.list()
    for i in range(d // l):
        if hlist[i] != 0:
            print('h not 0')
            return 1
    gamma, gamma_hash = get_gamma(t0, t1, t2, w)
    intt_factor = l * gamma
    tau = (intt_factor * t1 - gamma).mod(X**d + 1)
    vulp = (scalar(list((intt_factor * b[0][i] + b[1][i]) for i in range(BASELEN)), z) - c * (tau + t2 - h)).mod(X**d + 1)
    c_hash_prime = get_challenge_hash(gamma_hash, t3, vpp, h, vulp)
    if c_hash != c_hash_prime:
        print('c != c_prime')
        return 1
    return 0


def proof_correctness(v, Nc, public_seed):
    try_index = 0
    m_prime = INTT([1] * Nc + [0] * (l - Nc))
    B0, b = gen_public_b_with_extra(public_seed)
    seed = randombytes(SEEDLEN)
    nonce = 0
    g, nonce = random_poly_with_zeros(seed, nonce, d, d//l)
    m = INTT(v)
    r_seed = randombytes(SEEDLEN)
    t0, t1, r, _ = commit(B0, b[0], m, r_seed, 0)
    t2 = scalar(b[1], r) + g
    while True:
        print('iteration', try_index)
        try_index += 1
    
        y, nonce = uniform_poly(BASELEN, seed, nonce)
        w = matrix_vector(B0, y)

        gamma, gamma_hash = get_gamma(t0, t1, t2, w)
        t3 = (scalar(b[2], r) - (2 * m - m_prime) * scalar(b[0], y)).mod(X**d + 1)
        vpp = (scalar(b[2], y) + scalar(b[0], y)**2).mod(X**d + 1)
        intt_factor = gamma * l
        h = (g + intt_factor * m - gamma).mod(X**d + 1)
        vulp = scalar(list((intt_factor * b[0][i] + b[1][i]) for i in range(BASELEN)), y)

        c_hash = get_challenge_hash(gamma_hash, t3, vpp, h, vulp)
        c = get_challenge(c_hash)

        z = [0] * BASELEN
        for i in range(BASELEN):
            z[i] = (y[i] + c * r[i]).mod(X**d + 1)
        if not check_z_len(z):
            break
    return (h, c_hash, z), (t0, t1, t2, t3), r


def verify_correctness(proof, commitment, public_seed):
    m_prime = INTT([1] * Nc + [0] * (l - Nc))
    B0, b = gen_public_b_with_extra(public_seed)
    h, c_hash, z = proof
    c = get_challenge(c_hash)
    t0, t1, t2, t3 = commitment
    if check_z_len(z):
        print('check_z_len')
        return 1
    B0z = matrix_vector(B0, z)
    w = [0] * KAPPA
    for i in range(KAPPA):
        w[i] = (B0z[i] - c * t0[i]).mod(X**d + 1)
    f1 = (scalar(b[0], z) - c * t1).mod(X**d + 1)
    f2 = (scalar(b[0], z) - c * (t1 - m_prime)).mod(X**d + 1)
    f3 = (scalar(b[2], z) - c * t3).mod(X**d + 1)
    vpp = (f1 * f2 + f3).mod(X**d + 1)
    hlist = h.list()
    for i in range(d // l):
        if hlist[i] != 0:
            print('h not 0')
            return 1
    gamma, gamma_hash = get_gamma(t0, t1, t2, w)
    intt_factor = l * gamma
    tau = (intt_factor * t1 - gamma).mod(X**d + 1)
    vulp = (scalar(list((intt_factor * b[0][i] + b[1][i]) for i in range(BASELEN)), z) - c * (tau + t2 - h)).mod(X**d + 1)
    c_hash_prime = get_challenge_hash(gamma_hash, t3, vpp, h, vulp)
    if c_hash != c_hash_prime:
        print('c != c_prime')
        return 1
    return 0


if __name__ == '__main__':
    from params import l, Nc

    v = [0] * l
    v[4] = 1

    public_seed = b'-\xc2\xbd\xc1\x12\x94\xac\xd0f\xab~\x9f\x13\xb5\xac\xcaT\xbaFgD\xa6\x93\xd9\x92\xf2"\xb5\x006\x02\xa3'
    proof, commitment, r = proof_correctness(v, Nc, public_seed)
    ver_result = verify_correctness(proof, commitment, public_seed)
    if ver_result == 0:
        print('Verify is successfull')
    else:
        print('There is an error in verification')