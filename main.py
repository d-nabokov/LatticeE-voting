from params import Nv_max, Na
from voting_scheme import vote, tally_j, tally_all, verify
from voting_classes import BBoard
from utils import max_with_index


def emulate_voting(BB, votes, public_seed):
    print('creating votes')
    for i in range(Nv):
        vote(i, votes[i], public_seed, BB)
    print('creating tallies')
    for j in range(Na):
        tally_j(j, public_seed, BB)
    print('computing total results')
    res = tally_all(public_seed, BB)
    winner_votes, winner_index = max_with_index(res)
    print(f'results of voting:\n{res}\nWinner is candidate {winner_index} with {winner_votes} votes')
    
    ver_result = verify(res, public_seed, BB)
    if ver_result == 0:
        print('Voting is successfull')
    else:
        print('There is an error in voting')


if __name__ == '__main__':
    public_seed = b'\xc2\x84\x80y\xef\xfew\xaf\n\x03\x95h\xa1\xee\xda}D\xbf\x87\x10a\xf3\xc6\x92\xe7\xa3\xa3\x9dTR\tY'
    Nv = 6
    assert(Nv <= Nv_max)
    votes = [0, 1, 2, 1, 1, 1]
    CA = {}
    BB = BBoard(CA)
    emulate_voting(BB, votes, public_seed)
