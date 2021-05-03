from collections import defaultdict

class BBoard:
    def __init__(self, CA, authorities):
        self.authorities = authorities
        self.ballots = {}
        self.tallies = {}
        self.ballot_correctness = {}
        self.CA = CA
        
        
    def add_ballot(self, ballot):
        self.ballots[ballot.v_id] = ballot
        
        
    def get_ballot(self, v_id):
        return self.ballots[v_id]
    
    
    def all_ballots(self):
        for v_id, ballot in self.ballots.items():
            yield v_id, ballot


    def add_authority_j_ballot_correctness(self, a_id, ballot_correctness):
        self.ballot_correctness[a_id] = ballot_correctness


    def correct_ballots(self):
        assert(len(self.authorities) == len(self.ballot_correctness))
        is_correct = defaultdict(int)
        for a_id, ballot_correctness in self.ballot_correctness.items():
            for v_id, is_ok in ballot_correctness.items():
                is_correct[v_id] |= is_ok
        for v_id, ballot in self.ballots.items():
            if is_correct[v_id] == 0:
                yield v_id, ballot
            
            
    def add_authority_tally(self, tally):
        self.tallies[tally.a_id] = tally
        
    
    def all_tallies(self):
        for a_id, tally in self.tallies.items():
            yield a_id, tally
            
    

class Ballot:
    def __init__(self, v_id, vproof, com, additional_com, enc_r, signature):
        self.v_id = v_id
        self.vproof = vproof
        self.com = com
        self.additional_com = additional_com
        self.enc_r = enc_r
        self.signature = signature
        
        
class Tally:
    def __init__(self, a_id, amo_proof, amo_zero_proof, T, r, signature):
        self.a_id = a_id
        self.amo_proof = amo_proof
        self.amo_zero_proof = amo_zero_proof
        self.additional_com = T
        self.final_com_r = r
        self.signature = signature


class Voter:
    def __init__(self, v_id, v_pk):
        self.v_id = v_id
        self.v_pk = v_pk


class Authority:
    def __init__(self, a_id, a_pk):
        self.a_id = a_id
        self.a_pk = a_pk