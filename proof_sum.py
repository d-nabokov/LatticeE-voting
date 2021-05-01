from sage.all import *

from params import PP.kappa, PP.seedlen, PP.baselen, R
from public import gen_public_b
from proof_amo import proof_amo, proof_amo_to_zero, verify_amo, verify_amo_to_zero
from linear_alg import scalar
from utils import randombytes
from commit import commit


# S = (s_1, ..., s_p); T0 = (t0_1, ..., t0_p); T1 = (t1_1, ..., t1_p)
def sum_of_commitments(S, T0, T1, u, public_seed):
    p = len(S)
    max_layers = ceil(log(p, u))
    
    B0, b1 = gen_public_b(PP, public_seed)
    X = list(T1[i] - scalar(b1, S[i]) for i in range(p))
    
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
        for b in range(blocks):
            start = b * u
            end = min((b + 1) * u, max_index)
            m = sum(X[i] for i in range(start, end))
            X[b] = m
            t0, t1, r, nonce = commit(B0, b1, m, r_seed, nonce)
            
            t0_prime = list(t0[i] - sum(T0_amo[j + offset][i] for j in range(start, end)) for i in range(PP.kappa))
            t1_prime = t1 - sum(T1_amo[j + offset] for j in range(start, end))
            r_prime = list(r[i] - sum(S_amo[j + offset][i] for j in range(start, end)) for i in range(PP.baselen))
            
            T0_amo.append(t0)
            T1_amo.append(t1)
            S_amo.append(r)
            T0_amo_zero.append(t0_prime)
            T1_amo_zero.append(t1_prime)
            S_amo_zero.append(r_prime)
        offset += max_index
        max_index = blocks
    
    amo_proof = proof_amo(Matrix(R, S_amo).transpose(), Matrix(R, T0_amo).transpose(), len(T1_amo), public_seed)
    amo_zero_proof = proof_amo_to_zero(
        Matrix(R, S_amo_zero).transpose(), Matrix(R, T0_amo_zero).transpose(), vector(R, T1_amo_zero),
        len(T1_amo_zero), public_seed)
    return amo_proof, amo_zero_proof, (T0_amo[p:], T1_amo[p:]), S_amo[-1]


def verify_sum_of_commitments(amo_proof, amo_zero_proof, T_orig, T, u, public_seed):
    T0_orig, T1_orig = T_orig
    p = len(T1_orig)
    max_layers = ceil(log(p, u))
    T0, T1 = T
    T0_amo = T0_orig + T0
    T1_amo = T1_orig + T1
    if verify_amo(amo_proof, Matrix(R, T0_amo).transpose(), len(T1_amo), public_seed) != 0:
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
            t1_prime = t1 - sum(T1_amo[j + offset] for j in range(start, end))
            
            T0_amo_zero.append(t0_prime)
            T1_amo_zero.append(t1_prime)
        
        offset_blocks += blocks
        offset += max_index
        max_index = blocks
    if verify_amo_to_zero(amo_zero_proof, Matrix(R, T0_amo_zero).transpose(), vector(R, T1_amo_zero), \
                          len(T1_amo_zero), public_seed) != 0:
        print('Verification of amortized proof of opening to zero failed')
        return 1
    return 0


if __name__ == '__main__':
    public_seed = b'\xa3\xe4\xf3\xf9Gl\xb69\xe2\xff~\x02\x087I\x18\x9a\x08\x88\x15\xe1\x83\x02\x7fP\xd2\x13-\xa1\xb5.\x88'
    from params import u, Nv
    from proof_amo import gen_random_commitments
    S, T0, T1 = gen_random_commitments(Nv, public_seed)
    S, T0, T1 = list(S.transpose()), list(T0.transpose()), list(T1)
    amo_proof, amo_zero_proof, T, r = sum_of_commitments(S, T0, T1, u, public_seed)
    ver_result = verify_sum_of_commitments(amo_proof, amo_zero_proof, (T0, T1), T, u, public_seed)
    if ver_result == 0:
        print('Verify is successfull')
    else:
        print('There is an error in verification')