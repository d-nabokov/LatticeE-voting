from sage.all import *

from public import gen_public_b
from proof_amo import proof_amo, proof_amo_to_zero, verify_amo, verify_amo_to_zero
from linear_alg import scalar
from utils import randombytes
from commit import commit


# S = (s_1, ..., s_p); T0 = (t0_1, ..., t0_p); T1 = (t1_1, ..., t1_p)
def sum_of_commitments(PP, S, T0, T1, public_seed):
    d = PP.d
    X = PP.X
    u = PP.u
    npoly = PP.npoly

    p = len(S)
    max_layers = ceil(log(p, u))
    
    B0, b = gen_public_b(PP, public_seed)
    X_poly = list(list(T1[i][j] - scalar(b[j], S[i], X, d) for j in range(npoly)) for i in range(p))
    
    S_amo = []
    T0_amo = []
    T1_amo = []
    S_amo_zero = []
    T0_amo_zero = []
    T1_amo_zero = []
    for i in range(p):
        S_amo.append(S[i])
        T0_amo.append(T0[i])
        T1_amo.append(T1[i])
    r_seed = randombytes(PP.seedlen)
    nonce = 0
    max_index = p
    offset = 0
    for layer in range(1, max_layers+1):
        blocks = ceil(p / u**layer)
        for block in range(blocks):
            start = block * u
            end = min((block + 1) * u, max_index)
            m = list(sum(X_poly[i][j] for i in range(start, end)) for j in range(npoly))
            X_poly[block] = m
            t0, t1, r, nonce = commit(PP, B0, b, m, r_seed, nonce)
            
            t0_prime = list(t0[i] - sum(T0_amo[j + offset][i] for j in range(start, end)) for i in range(PP.kappa))
            t1_prime = list(t1[i] - sum(T1_amo[j + offset][i] for j in range(start, end)) for i in range(npoly))
            r_prime = list(r[i] - sum(S_amo[j + offset][i] for j in range(start, end)) for i in range(PP.baselen))
            
            T0_amo.append(t0)
            T1_amo.append(t1)
            S_amo.append(r)
            T0_amo_zero.append(t0_prime)
            T1_amo_zero.append(t1_prime)
            S_amo_zero.append(r_prime)
        offset += max_index
        max_index = blocks
    
    R = PP.R
    amo_proof = proof_amo(PP, Matrix(R, S_amo).transpose(), Matrix(R, T0_amo).transpose(), len(T1_amo), public_seed)
    amo_zero_proof = proof_amo_to_zero(
        PP, Matrix(R, S_amo_zero).transpose(), Matrix(R, T0_amo_zero).transpose(), Matrix(R, T1_amo_zero).transpose(),
        len(T1_amo_zero), public_seed)
    return amo_proof, amo_zero_proof, (T0_amo[p:], T1_amo[p:]), S_amo[-1]


def verify_sum_of_commitments(PP, amo_proof, amo_zero_proof, T_orig, T, public_seed):
    d = PP.d
    X = PP.X
    u = PP.u
    R = PP.R

    T0_orig, T1_orig = T_orig
    p = len(T1_orig)
    max_layers = ceil(log(p, u))
    T0, T1 = T
    T0_amo = T0_orig + T0
    T1_amo = T1_orig + T1
    if verify_amo(PP, amo_proof, Matrix(R, T0_amo).transpose(), len(T1_amo), public_seed) != 0:
        print('Verification of amortized proof failed')
        return 1
    
    T0_amo_zero = []
    T1_amo_zero = []
    max_index = p
    offset = 0
    offset_blocks = 0
    for layer in range(1, max_layers+1):
        blocks = ceil(p / u**layer)
        for b in range(blocks):
            start = b * u
            end = min((b + 1) * u, max_index)
            
            t0 = T0[offset_blocks + b]
            t1 = T1[offset_blocks + b]
            
            t0_prime = list(t0[i] - sum(T0_amo[j + offset][i] for j in range(start, end)) for i in range(PP.kappa))
            t1_prime = list(t1[i] - sum(T1_amo[j + offset][i] for j in range(start, end)) for i in range(PP.npoly))
            
            T0_amo_zero.append(t0_prime)
            T1_amo_zero.append(t1_prime)
        
        offset_blocks += blocks
        offset += max_index
        max_index = blocks
    if verify_amo_to_zero(PP, amo_zero_proof, Matrix(R, T0_amo_zero).transpose(), Matrix(R, T1_amo_zero).transpose(), \
                          len(T1_amo_zero), public_seed) != 0:
        print('Verification of amortized proof of opening to zero failed')
        return 1
    return 0


if __name__ == '__main__':
    public_seed = b'\xa3\xe4\xf3\xf9Gl\xb69\xe2\xff~\x02\x087I\x18\x9a\x08\x88\x15\xe1\x83\x02\x7fP\xd2\x13-\xa1\xb5.\x88'
    from params import PublicParams
    from proof_amo import gen_random_commitments
    PP = PublicParams(2, 5, 10)
    Nv = 4
    S, T0, T1 = gen_random_commitments(PP, Nv, public_seed)
    S, T0, T1 = list(S.transpose()), list(T0.transpose()), list(T1.transpose())
    amo_proof, amo_zero_proof, T, r = sum_of_commitments(PP, S, T0, T1, public_seed)
    ver_result = verify_sum_of_commitments(PP, amo_proof, amo_zero_proof, (T0, T1), T, public_seed)
    if ver_result == 0:
        print('Verify is successfull')
    else:
        print('There is an error in verification')