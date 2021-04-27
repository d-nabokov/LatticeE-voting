from Crypto.Hash import SHAKE128
from sage.stats.distributions.discrete_gaussian_integer import DiscreteGaussianDistributionIntegerSampler
from sage.all import *

from params import d, LOGDELTA1, DELTA1, q, q_bits, P, R, X, sigma


def random_poly(seed, nonce, degree):
    return random_poly_with_zeros(seed, nonce, degree, 0)


def random_poly_with_zeros(seed, nonce, degree, zeros):
    shake = SHAKE128.new()
    shake.update(int(nonce).to_bytes(8, byteorder='big'))
    shake.update(seed)
    rnd = shake.read(int(4 * degree))
    rnd_cnt = 0
    lst = [0] * degree
    i = zeros
    while i < degree:
        if rnd_cnt < degree:
            r = int.from_bytes(rnd[rnd_cnt*4 : (rnd_cnt + 1)*4], byteorder='big')
            rnd_cnt += 1
        else:
            r = int.from_bytes(shake.read(int(4)), byteorder='big')
        r = r & ((1 << q_bits) - 1)
        if r >= q:
            continue
        else:
            lst[i] = r
            i += 1
    return P(lst), nonce + 1


def _chi_map(c):
    return c - (3 & -(c >> 1))


def _chi_poly_single(seed, nonce):
    shake = SHAKE128.new()
    shake.update(int(nonce).to_bytes(8, byteorder='big'))
    shake.update(seed)
    samples_len = d
    rnd = shake.read(int(samples_len * 4 // 8))
    
    arr = [0]*(samples_len)
    a = [0] * 4
    for i in range(samples_len // 2):
        r = rnd[i]
        for j in range(2):
            for k in range(4):
                a[k] = r & 1
                r >>= 1
            c = (a[0] + a[1] - a[2] - a[3]) % 3
            arr[i * 2 + j] = _chi_map(c)
    return P(arr), nonce + 1

    
def chi_poly(n, seed, nonce):
    res = []
    for i in range(n):
        p, nonce = _chi_poly_single(seed, nonce)
        res.append(p)
    return res, nonce


def challenge(seed):
    shake = SHAKE128.new()
    shake.update(seed)
    rnd = shake.read(int(d * 2 // 8))
    
    arr = [0] * d
    a = [0] * 2
    for i in range(d // 4):
        r = rnd[i]
        for j in range(4):
            for k in range(2):
                a[k] = r & 1
                r >>= 1
            c = a[0] - a[1]
            arr[i * 4 + j] = c
    return P(arr)


def _uniform_poly_single(seed, nonce):
    shake = SHAKE128.new()
    shake.update(int(nonce).to_bytes(8, byteorder='big'))
    shake.update(seed)
    # absolutely not effective
    deltabytes = (LOGDELTA1 + 1 + 7) // 8
    assert(LOGDELTA1 + 1 <= deltabytes * 8)
    rnd = shake.read(int(deltabytes * d))

    samples_len = d
    arr = [0]*(samples_len)
    i = 0
    rnd_offset = 0
    while i < samples_len:
        if rnd_offset < len(rnd):
            r = int.from_bytes(rnd[rnd_offset : rnd_offset + deltabytes], byteorder='big')
            rnd_offset += deltabytes
        else:
            r = int.from_bytes(shake.read(int(deltabytes)), byteorder='big')
        r = (r & ((DELTA1 << 1) - 1)) - DELTA1
        if r != -DELTA1:
            arr[i] = r
            i += 1
    return P(arr), nonce + 1


def uniform_poly(n, seed, nonce):
    res = []
    for i in range(n):
        p, nonce = _uniform_poly_single(seed, nonce)
        res.append(p)
    return res, nonce


def challenge_amo(seed, p, N):
    shake = SHAKE128.new()
    shake.update(seed)
    # absolutely not effective
    val_bytes = 2
    rnd = shake.read(int(val_bytes * 2 * p * N))
    rnd_offset = 0
    
    c = []
    for i in range(p):
        row = [0] * N
        for j in range(N):
            while True:
                if rnd_offset < len(rnd):
                    r = int.from_bytes(rnd[rnd_offset : rnd_offset + val_bytes], byteorder='big')
                    rnd_offset += val_bytes
                else:
                    r = int.from_bytes(shake.read(int(val_bytes)), byteorder='big')
                r & 511
                if r <= 257:
                    if r == 0:
                        row[j] = 0
                    else:
                        row[j] = (-1)**(r & 1) * X**((r - 1) // 2)
                    break
        c.append(row)
    return c


def discrete_gaussian_y(m, n):
    return Matrix(R, m, n, lambda i, j: R(list(D() for _ in range(d))))


D = DiscreteGaussianDistributionIntegerSampler(sigma=sigma)