from Crypto.Hash import SHAKE128
from sage.all import *

from params import PP.kappa, PP.seedlen, PP.baselen, d, R, p, N, M, BETA2, sigma
from public import gen_public_b
from utils import poly_to_bytes, rejection_sampling, randombytes
from linear_alg import matrix_vector, scalar
from random_polynomials import challenge_amo, random_poly, discrete_gaussian_y, chi_poly
from commit import commit
from ring import signed_zq


def get_challenge_hash_amo(T, W, p):
    shake = SHAKE128.new()
    for i in range(PP.kappa):
        for j in range(p):
            shake.update(poly_to_bytes(T[i][j]))    
    for i in range(PP.kappa):
        for j in range(N):
            shake.update(poly_to_bytes(W[i][j]))
    c_hash = shake.read(int(PP.seedlen))
    return c_hash


def get_challenge_amo(c_hash, p):
    return challenge_amo(c_hash, p, N)


# def check_Z_len_infty(Z):
#     border = DELTA1 - BETA1
#     norm = inf_norm_matr(Z)
#     if norm >= border:
#         return 1
#     return 0


# def proof_amo_infty(S, T, public_seed):
#     B0, b1 = gen_public_b(PP, public_seed)
#     try_index = 0
#     seed = randombytes(PP.seedlen)
#     nonce = 0
#     while True:
#         print('iteration', try_index)
#         try_index += 1
    
#         Y = []
#         for i in range(N):
#             Y.append(uniform_poly(PP.baselen, seed, nonce))
#             nonce += PP.baselen
#         Y = Matrix(R, Y).transpose()
#         W = Matrix(R, B0) * Y
        
#         C, c_hash = get_challenge_amo(T, W)
#         C = Matrix(R, C)
#         Z = Y + S * C
#         print(f'||SC|| = {inf_norm_matr(S * C)}')
#         print(f'||Y|| = {inf_norm_matr(Y)}')
#         print(f'||Z|| = {inf_norm_matr(Z)}')
#         if not check_Z_len(Z):
#             break
            
#     return (C, Z)
        
        
# def verify_amo_infty(proof, T, public_seed):
#     C, Z = proof
#     B0, b1 = gen_public_b(PP, public_seed)
    
#     if check_Z_len(Z):
#         print('check_Z_len')
#         return 1
    
#     W = Matrix(R, B0) * Z - T * C
#     C_prime, c_hash = get_challenge_amo(T, W)
#     if C != Matrix(R, C_prime):
#         print('C != C_prime')
#         return 1
#     return 0


def check_Z_len(Z):
    norm = [0] * N
    for row in Z:
        for j, val in enumerate(row):
            for coef in val:
                norm[j] += signed_zq(coef)**2
    
    for i in range(N):
        norm[i] = sqrt(norm[i])
        if norm[i] > BETA2:
            return 1
    return 0


def gen_randomness_matrix(p):
    r_seed = randombytes(PP.seedlen)
    S = []
    for i in range(p):
        r = chi_poly(PP.baselen, r_seed, i)
        S.append(r)
    return Matrix(R, S).transpose()


def proof_amo(S, T, p, public_seed):
    B0, b1 = gen_public_b(PP, public_seed)
    try_index = 0
#     seed = randombytes(PP.seedlen)
#     nonce = 0
    while True:
        print('iteration', try_index)
        try_index += 1
    
        Y = discrete_gaussian_y(PP.baselen, N)
        W = Matrix(R, B0) * Y
        
        c_hash = get_challenge_hash_amo(T, W, p)
        C = get_challenge_amo(c_hash, p)
        C = Matrix(R, C)
        SC = S * C
        Z = Y + SC
        if rejection_sampling(Z, SC, sigma, M):
            break
            
    return (c_hash, Z)
        
        
def verify_amo(proof, T, p, public_seed):
    c_hash, Z = proof
    C = get_challenge_amo(c_hash, p)
    C = Matrix(R, C)
    B0, b1 = gen_public_b(PP, public_seed)
    
    if check_Z_len(Z):
        print('check_Z_len')
        return 1
    
    W = Matrix(R, B0) * Z - T * C
    c_hash_prime = get_challenge_hash_amo(T, W, p)
    if c_hash != c_hash_prime:
        print('C != C_prime')
        return 1
    return 0


def get_challenge_hash_amo_to_zero(T0, T1, W0, W1, p):
    shake = SHAKE128.new()
    for i in range(PP.kappa):
        for j in range(p):
            shake.update(poly_to_bytes(T0[i][j]))
    for i in range(p):
        shake.update(poly_to_bytes(T1[i]))
    for i in range(PP.kappa):
        for j in range(N):
            shake.update(poly_to_bytes(W0[i][j]))
    for i in range(N):
        shake.update(poly_to_bytes(W1[i]))
    c_hash = shake.read(int(PP.seedlen))
    return c_hash


def get_challenge_amo_to_zero(c_hash, p):
    return challenge_amo(c_hash, p, N)


def proof_amo_to_zero(S, T0, T1, p, public_seed):
    B0, b1 = gen_public_b(PP, public_seed)
    try_index = 0
    while True:
        print('iteration', try_index)
        try_index += 1
    
        Y = discrete_gaussian_y(PP.baselen, N)
        W0 = Matrix(R, B0) * Y
        W1 = vector(R, b1) * Y
        
        c_hash = get_challenge_hash_amo_to_zero(T0, T1, W0, W1, p)
        C = get_challenge_amo_to_zero(c_hash, p)
        C = Matrix(R, C)
        SC = S * C
        Z = Y + SC
        if rejection_sampling(Z, SC, sigma, M):
            break
   
    return (c_hash, Z)
        
        
def verify_amo_to_zero(proof, T0, T1, p, public_seed):
    c_hash, Z = proof
    C = get_challenge_amo_to_zero(c_hash, p)
    C = Matrix(R, C)
    B0, b1 = gen_public_b(PP, public_seed)
    
    if check_Z_len(Z):
        print('check_Z_len')
        return 1
    
    W0 = Matrix(R, B0) * Z - T0 * C
    W1 = vector(R, b1) * Z - T1 * C
    
    c_hash_prime = get_challenge_hash_amo_to_zero(T0, T1, W0, W1, p)
    if c_hash != c_hash_prime:
        print('C != C_prime')
        return 1
    return 0


def gen_random_commitments(p, public_seed):
    B0, b1 = gen_public_b(PP, public_seed)
    
    r_seed = randombytes(PP.seedlen)
    seed = randombytes(PP.seedlen)
    S = []
    T0 = []
    T1 = []
    nonce = 0
    for i in range(p):
        xi, nonce = random_poly(seed, nonce, d)
        t0, t1, r, nonce = commit(B0, b1, xi, r_seed, nonce)
        S.append(r)
        T0.append(t0)
        T1.append(t1)
    return Matrix(R, S).transpose(), Matrix(R, T0).transpose(), vector(R, T1)


if __name__ == '__main__':
    print(randombytes(PP.seedlen))
    public_seed = b'\xf3\xe0\xf0\n\x17\x02\xd3\xee\xd3\xbd{D\xff\x19\xf5b\x98\xca\xdf\xc0M\xe8\x12\xbe\xc3\xc4a1\xd6\xe1\xf2\xba'

    S, T0, T1 = gen_random_commitments(p, public_seed)
    proof = proof_amo(S, T0, p, public_seed)
    ver_result = verify_amo(proof, T0, p, public_seed)
    if ver_result == 0:
        print('Verify is successfull')
    else:
        print('There is an error in verification')