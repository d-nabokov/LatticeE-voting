from Crypto.Hash import SHAKE128
from sage.all import *

from public import gen_public_b
from utils import poly_to_bytes, rejection_sampling_matrix, randombytes
from linear_alg import matrix_vector, scalar, inf_norm_matr
from random_polynomials import challenge_amo, random_poly, discrete_gaussian_y, chi_poly
from commit import commit
from ring import signed_zq


def get_challenge_hash_amo(PP, T, W, p):
    shake = SHAKE128.new()
    for i in range(PP.kappa):
        for j in range(p):
            shake.update(poly_to_bytes(T[i][j]))    
    for i in range(PP.kappa):
        for j in range(PP.amo_n):
            shake.update(poly_to_bytes(W[i][j]))
    c_hash = shake.read(int(PP.seedlen))
    return c_hash


def get_challenge_amo(PP, c_hash, p):
    return challenge_amo(PP, c_hash, p)


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


def check_Z_len(PP, Z):
    norm = [0] * PP.amo_n
    for row in Z:
        for j, val in enumerate(row):
            for coef in val:
                norm[j] += signed_zq(coef, PP.q)**2
    
    for i in range(PP.amo_n):
        norm[i] = sqrt(norm[i])
        if norm[i] > PP.beta2:
            return 1
    return 0


def gen_randomness_matrix(PP, p):
    R = PP.R

    r_seed = randombytes(PP.seedlen)
    S = []
    for i in range(p):
        r = chi_poly(PP.baselen, r_seed, i)
        S.append(r)
    return Matrix(R, S).transpose()


def proof_amo(PP, S, T, p, public_seed):
    R = PP.R

    B0, b = gen_public_b(PP, public_seed)
    while True:
        Y = discrete_gaussian_y(PP, PP.baselen, PP.amo_n, PP.sigma2)
        W = Matrix(R, B0) * Y
        
        c_hash = get_challenge_hash_amo(PP, T, W, p)
        C = get_challenge_amo(PP, c_hash, p)
        C = Matrix(R, C)
        SC = S * C
        Z = Y + SC
        if rejection_sampling_matrix(Z, SC, PP.sigma2, PP.average_rejection_tries2, PP.q) \
            and inf_norm_matr(Z, PP.q) < PP.inf_bound2:
            break
            
    return (c_hash, Z)
        
        
def verify_amo(PP, proof, T, p, public_seed):
    R = PP.R

    c_hash, Z = proof
    C = get_challenge_amo(PP, c_hash, p)
    C = Matrix(R, C)
    B0, b = gen_public_b(PP, public_seed)
    
    if check_Z_len(PP, Z):
        return 1
    
    W = Matrix(R, B0) * Z - T * C
    c_hash_prime = get_challenge_hash_amo(PP, T, W, p)
    if c_hash != c_hash_prime:
        return 1
    return 0


def get_challenge_hash_amo_to_zero(PP, T0, T1, W0, W1, p):
    shake = SHAKE128.new()
    for i in range(PP.kappa):
        for j in range(p):
            shake.update(poly_to_bytes(T0[i][j]))
    for i in range(PP.npoly):
        for j in range(p):
            shake.update(poly_to_bytes(T1[i][j]))
    for i in range(PP.kappa):
        for j in range(PP.amo_n):
            shake.update(poly_to_bytes(W0[i][j]))
    for i in range(PP.npoly):
        for j in range(PP.amo_n):
            shake.update(poly_to_bytes(W1[i][j]))
    c_hash = shake.read(int(PP.seedlen))
    return c_hash


def get_challenge_amo_to_zero(PP, c_hash, p):
    return challenge_amo(PP, c_hash, p)


def proof_amo_to_zero(PP, S, T0, T1, p, public_seed):
    R = PP.R

    B0, b = gen_public_b(PP, public_seed)
    while True:
        Y = discrete_gaussian_y(PP, PP.baselen, PP.amo_n, PP.sigma2)
        W0 = Matrix(R, B0) * Y
        W1 = Matrix(R, list(vector(R, b[i]) * Y for i in range(PP.npoly)))
        
        c_hash = get_challenge_hash_amo_to_zero(PP, T0, T1, W0, W1, p)
        C = get_challenge_amo_to_zero(PP, c_hash, p)
        C = Matrix(R, C)
        SC = S * C
        Z = Y + SC
        if rejection_sampling_matrix(Z, SC, PP.sigma2, PP.average_rejection_tries2, PP.q) \
            and inf_norm_matr(Z, PP.q) < PP.inf_bound2:
            break
   
    return (c_hash, Z)
        
        
def verify_amo_to_zero(PP, proof, T0, T1, p, public_seed):
    R = PP.R

    c_hash, Z = proof
    C = get_challenge_amo_to_zero(PP, c_hash, p)
    C = Matrix(R, C)
    B0, b = gen_public_b(PP, public_seed)
    
    if check_Z_len(PP, Z):
        return 1
    
    W0 = Matrix(R, B0) * Z - T0 * C
    W1 = Matrix(R, list(vector(R, b[i]) * Z - T1[i] * C for i in range(PP.npoly)))
    
    c_hash_prime = get_challenge_hash_amo_to_zero(PP, T0, T1, W0, W1, p)
    if c_hash != c_hash_prime:
        return 1
    return 0


def _gen_random_commitments_impl(PP, p, public_seed, m_func):
    R = PP.R

    B0, b = gen_public_b(PP, public_seed)
    
    r_seed = randombytes(PP.seedlen)
    seed = randombytes(PP.seedlen)
    S = []
    T0 = []
    T1 = []
    nonce = 0
    for i in range(p):
        xi = []
        for i in range(PP.npoly):
            x, nonce = m_func(PP, seed, nonce, PP.d)
            xi.append(x)
        t0, t1, r, nonce = commit(PP, B0, b, xi, r_seed, nonce)
        S.append(r)
        T0.append(t0)
        T1.append(t1)
    return Matrix(R, S).transpose(), Matrix(R, T0).transpose(), Matrix(R, T1).transpose()


def gen_random_commitments(PP, p, public_seed):
    return _gen_random_commitments_impl(PP, p, public_seed, random_poly)


def gen_random_commitments_to_zero(PP, p, public_seed):
    return _gen_random_commitments_impl(PP, p, public_seed, lambda a, b, c, d: (0, 0))


if __name__ == '__main__':
    from params import PublicParams
    public_seed = b'\xf3\xe0\xf0\n\x17\x02\xd3\xee\xd3\xbd{D\xff\x19\xf5b\x98\xca\xdf\xc0M\xe8\x12\xbe\xc3\xc4a1\xd6\xe1\xf2\xba'

    PP = PublicParams(2, 5, 10)
    Nv = 7
    p = PP.number_of_authority_commitments(Nv)
    S, T0, T1 = gen_random_commitments(PP, p, public_seed)
    proof = proof_amo(PP, S, T0, p, public_seed)
    ver_result = verify_amo(PP, proof, T0, p, public_seed)
    if ver_result == 0:
        print('Verify amo is successfull')
    else:
        print('There is an error in verification amo')
    S, T0, T1 = gen_random_commitments_to_zero(PP, p - Nv, public_seed)
    amo_zero_proof = proof_amo_to_zero(PP, S, T0, T1, p - Nv, public_seed)
    ver_result_zero = verify_amo_to_zero(PP, amo_zero_proof, T0, T1, p - Nv, public_seed)
    if ver_result_zero == 0:
        print('Verify amo to zero is successfull')
    else:
        print('There is an error in verification amo to zero')
