from random_polynomials import random_poly
from commit import commit, commit_with_r
from utils import randombytes, m_from_vote_arr
from linear_alg import inf_norm_vect, scalar, matrix_vector
from public import gen_public_b
from ring import NTT, INTT
from proof_v import proof_v, verify_v
from proof_sum import sum_of_commitments, verify_sum_of_commitments
from voting_classes import Tally, Ballot


def _secret_share_value(PP, val, n, seed, nonce):
    res = []
    for i in range(n):
        res.append([0] * PP.npoly)
    for j in range(PP.npoly):
        x = [0] * n
        for i in range(n - 1):
            x[i], nonce = random_poly(PP, seed, nonce, PP.d)
        x[n - 1] = val[j] - sum(x[i] for i in range(n - 1))
        for i in range(n):
            res[i][j] = x[i]
    return res, nonce


def vote(PP, v_id, vote_arr, public_seed, BB):
    seed = randombytes(PP.seedlen)
    r_seed = randombytes(PP.seedlen)
    nonce = 0
    
    B0, b = gen_public_b(PP, public_seed)
    m = m_from_vote_arr(PP, vote_arr)
    x, _ = _secret_share_value(PP, m, PP.Na, seed, 0)
    
    S = []
    T0 = []
    T1 = []
    for i in range(PP.Na):
        t0, t1, r, nonce = commit(PP, B0, b, x[i], r_seed, nonce)
        S.append(r)
        T0.append(t0)
        T1.append(t1)
    
    r = list(sum(S[j][i] for j in range(PP.Na)) for i in range(PP.baselen))
    t0, t1 = commit_with_r(PP, B0, b, m, r)
    vproof, additional_com = proof_v(PP, t0, t1, r, m, public_seed)
    
    # We should sign ballots
    signature = None
    
    # S should be encrypted
    ballot = Ballot(v_id, vproof, (T0, T1), additional_com, S, signature)
    BB.add_ballot(ballot)
    

def _check_ballot_sig(ballot, CA):
    return 1


def _check_tally_sig(tally, CA):
    return 1


def extract_from_ballot_by_j(PP, v_id, a_id, BB):
    ballot = BB.get_ballot(v_id)
    # Should check signature
    ballot.signature
    # Should decrypt using authority's secret key
    r = ballot.enc_r[a_id]
    if inf_norm_vect(r, PP.q) > PP.beta_commit_infty:
        print('r has invalid values')
        return 1
    T0, T1 = ballot.com
    return r, T0[a_id], T1[a_id]


def tally_j(PP, a_id, public_seed, BB):
    S = []
    T0 = []
    T1 = []
    for v_id, ballot in BB.all_ballots():
        r, t0, t1 = extract_from_ballot_by_j(PP, v_id, a_id, BB)
        S.append(r)
        T0.append(t0)
        T1.append(t1)
    amo_proof, amo_zero_proof, T, r = sum_of_commitments(PP, S, T0, T1, public_seed)
    # We should sign tallies
    signature = None
    BB.add_authority_tally(Tally(a_id, amo_proof, amo_zero_proof, T, r, signature))
    
    
def _extract_authority_result(PP, B0, b, tally):
    r = tally.final_com_r
    T0, T1 = tally.additional_com
    t0 = T0[-1]
    t1 = T1[-1]
    B0r = matrix_vector(B0, r, PP.X, PP.d)
    zero = list(t0[i] - B0r[i] for i in range(PP.kappa))
    x = list((t1[i] - scalar(b[i], r, PP.X, PP.d)) for i in range(PP.npoly))
    return zero, x


def _check_zero(PP, zero):
    for i in range(PP.kappa):
        if zero[i] != 0:
            return 1
    return 0


def _res_to_list_of_candidates(PP, res):
    l = list(val for x in res for val in NTT(PP, x))
    return l[:PP.Nc]
    
    
def tally_all(PP, public_seed, BB):
    B0, b = gen_public_b(PP, public_seed)
    res = [0] * PP.npoly
    for a_id, tally in BB.all_tallies():
        zero, x = _extract_authority_result(PP, B0, b, tally)
        if _check_zero(PP, zero):
            raise Exception(f'autority {a_id} is cheating!')
        
        for i in range(PP.npoly):
            res[i] += x[i]
    return _res_to_list_of_candidates(PP, res)


def verify(PP, result, public_seed, BB):
    B0, b = gen_public_b(PP, public_seed)
    T0_a = []
    T1_a = []
    for i in range(PP.Na):
        T0_a.append(list())
        T1_a.append(list())
    
    for v_id, ballot in BB.all_ballots():
        T0, T1 = ballot.com
        for i in range(PP.Na):
            T0_a[i].append(T0[i])
            T1_a[i].append(T1[i])
        t0 = list(sum(T0[j][i] for j in range(PP.Na)) for i in range(PP.kappa))
        t1 = list(sum(T1[j][i] for j in range(PP.Na)) for i in range(PP.npoly))
        if verify_v(PP, ballot.vproof, (t0, t1), ballot.additional_com, public_seed):
            print(f'voter {v_id} is cheating!')
            return 1
    j = 0
    res = [0] * PP.npoly
    for a_id, tally in BB.all_tallies():
        if verify_sum_of_commitments(PP, tally.amo_proof, tally.amo_zero_proof,\
                                        (T0_a[j], T1_a[j]), tally.additional_com, public_seed):
            print(f'autority {a_id} is cheating!')
            return 1
        j += 1
        zero, x = _extract_authority_result(PP, B0, b, tally)
        if _check_zero(PP, zero):
            print(f'autority {a_id} is cheating!')
            return 1
        for i in range(PP.npoly):
            res[i] += x[i]
    if _res_to_list_of_candidates(PP, res) != result:
        print('results are wrong!')
        return 1
    return 0