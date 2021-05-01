from Crypto.Hash import SHAKE128

from ring import INTT
from public import gen_public_b_with_extra
from utils import poly_to_bytes, randombytes
from linear_alg import inf_norm, matrix_vector, scalar
from random_polynomials import challenge, random_poly, random_poly_with_zeros, uniform_poly


def get_gamma(PP, t0, t1, t2, w):
    shake = SHAKE128.new()
    for i in range(PP.kappa):
        shake.update(poly_to_bytes(t0[i]))
    shake.update(poly_to_bytes(t1))
    shake.update(poly_to_bytes(t2))
    for i in range(PP.kappa):
        shake.update(poly_to_bytes(w[i]))
    gamma_hash = shake.read(int(PP.seedlen))
    gamma, _ = random_poly(PP, gamma_hash, 0, PP.d//PP.l)
    return gamma, gamma_hash


def get_challenge_hash(PP, gamma_hash, t3, vpp, h, vulp):
    shake = SHAKE128.new()
    shake.update(gamma_hash)
    shake.update(poly_to_bytes(t3))
    shake.update(poly_to_bytes(vpp))
    shake.update(poly_to_bytes(h))
    shake.update(poly_to_bytes(vulp))
    c_hash = shake.read(int(PP.seedlen))
    return c_hash


def get_challenge(PP, c_hash):
    return challenge(PP, c_hash)


def check_z_len(PP, z):
    border = PP.delta1 - PP.beta1
    for i in range(PP.baselen):
        if inf_norm(z[i], PP.q) >= border:
            return 1
    return 0


def proof_v(PP, t0, t1, r, m, public_seed):
    d = PP.d
    l = PP.l
    X = PP.X

    # try_index = 0
    m_prime = INTT(PP, [1] * PP.Nc + [0] * (l - PP.Nc))
    B0, b = gen_public_b_with_extra(PP, public_seed)
    seed = randombytes(PP.seedlen)
    nonce = 0
    g, nonce = random_poly_with_zeros(PP, seed, nonce, d, d//l)
    t2 = scalar(b[1], r, X, d) + g
    while True:
        # print('iteration', try_index)
        # try_index += 1
    
        y, nonce = uniform_poly(PP, PP.baselen, seed, nonce)
        w = matrix_vector(B0, y, X, d)

        gamma, gamma_hash = get_gamma(PP, t0, t1, t2, w)
        t3 = (scalar(b[2], r, X, d) - (2 * m - m_prime) * scalar(b[0], y, X, d)).mod(X**d + 1)
        vpp = (scalar(b[2], y, X, d) + scalar(b[0], y, X, d)**2).mod(X**d + 1)
        intt_factor = gamma * l
        h = (g + intt_factor * m - gamma).mod(X**d + 1)
        vulp = scalar(list((intt_factor * b[0][i] + b[1][i]) for i in range(PP.baselen)), y, X, d)

        c_hash = get_challenge_hash(PP, gamma_hash, t3, vpp, h, vulp)
        c = get_challenge(PP, c_hash)

        z = [0] * PP.baselen
        for i in range(PP.baselen):
            z[i] = (y[i] + c * r[i]).mod(X**d + 1)
        if not check_z_len(PP, z):
            break
    return (h, c_hash, z), (t2, t3)


def verify_v(PP, proof, commitment, additional_com, public_seed):
    d = PP.d
    l = PP.l
    X = PP.X

    m_prime = INTT(PP, [1] * PP.Nc + [0] * (l - PP.Nc))
    B0, b = gen_public_b_with_extra(PP, public_seed)
    h, c_hash, z = proof
    c = get_challenge(PP, c_hash)
    t0, t1 = commitment
    t2, t3 = additional_com
    if check_z_len(PP, z):
        print('check_z_len')
        return 1
    B0z = matrix_vector(B0, z, X, d)
    w = [0] * PP.kappa
    for i in range(PP.kappa):
        w[i] = (B0z[i] - c * t0[i]).mod(X**d + 1)
    f1 = (scalar(b[0], z, X, d) - c * t1).mod(X**d + 1)
    f2 = (scalar(b[0], z, X, d) - c * (t1 - m_prime)).mod(X**d + 1)
    f3 = (scalar(b[2], z, X, d) - c * t3).mod(X**d + 1)
    vpp = (f1 * f2 + f3).mod(X**d + 1)
    hlist = h.list()
    for i in range(d // l):
        if hlist[i] != 0:
            print('h not 0')
            return 1
    gamma, gamma_hash = get_gamma(PP, t0, t1, t2, w)
    intt_factor = l * gamma
    tau = (intt_factor * t1 - gamma).mod(X**d + 1)
    vulp = (scalar(list((intt_factor * b[0][i] + b[1][i]) for i in range(PP.baselen)), z, X, d) - c * (tau + t2 - h)).mod(X**d + 1)
    c_hash_prime = get_challenge_hash(PP, gamma_hash, t3, vpp, h, vulp)
    if c_hash != c_hash_prime:
        print('c != c_prime')
        return 1
    return 0


if __name__ == '__main__':
    from commit import commit
    from params import PublicParams
    from public import gen_public_b
    public_seed = b'-\xc2\xbd\xc1\x12\x94\xac\xd0f\xab~\x9f\x13\xb5\xac\xcaT\xbaFgD\xa6\x93\xd9\x92\xf2"\xb5\x006\x02\xa3'

    PP = PublicParams(2, 5, 10)
    v = [0] * PP.l
    v[1] = 1
    m = INTT(PP, v)
    B0, b1 = gen_public_b(PP, public_seed)
    r_seed = randombytes(PP.seedlen)
    t0, t1, r, _ = commit(PP, B0, b1, m, r_seed, 0)

    proof, additional_com = proof_v(PP, t0, t1, r, m, public_seed)
    ver_result = verify_v(PP, proof, (t0, t1), additional_com, public_seed)
    if ver_result == 0:
        print('Verify is successfull')
    else:
        print('There is an error in verification')

    print('Trying negative scenarios')
    for v in ([1]*2 + [0]*(PP.l - 2), [2] + [0]*(PP.l - 1)):
        m = INTT(PP, v)
        B0, b1 = gen_public_b(PP, public_seed)
        r_seed = randombytes(PP.seedlen)
        t0, t1, r, _ = commit(PP, B0, b1, m, r_seed, 0)

        proof, additional_com = proof_v(PP, t0, t1, r, m, public_seed)
        ver_result = verify_v(PP, proof, (t0, t1), additional_com, public_seed)
        assert(ver_result == 1)
