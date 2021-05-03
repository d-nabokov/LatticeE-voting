from params import PublicParams
from voting_scheme import vote, testballots, tally_j, tally_all, verify
from voting_classes import BBoard
from utils import max_with_index, reset_all


def emulate_voting(PP, BB, votes, public_seed):
    Nv = len(votes)
    print('creating votes')
    for i in range(Nv):
        print(f'\tvoter {i} casting {votes[i]}')
        vote(PP, i, votes[i], public_seed, BB)
    print('testing ballots')
    for j in range(PP.Na):
        testballots(PP, j, public_seed, BB)
    print('creating tallies')
    for j in range(PP.Na):
        print(f'\tauthority {j} creating tally')
        tally_j(PP, j, public_seed, BB)
    print('computing total results')
    res = tally_all(PP, public_seed, BB)
    winner_votes, winner_indexes = max_with_index(res)
    winner_indexes = ', '.join(str(x) for x in winner_indexes)
    print(f'results of voting:\n{res}\nWinner is candidate(s) "{winner_indexes}" with {winner_votes} votes')
    
    ver_result = verify(PP, res, public_seed, BB)
    if ver_result == 0:
        print('Voting is successfull')
    else:
        print('There is an error in voting')


def vote_from_candidate(Nc, candidate):
    v = [0] * Nc
    v[candidate] = 1
    return v


if __name__ == '__main__':
    public_seed = b'\xc2\x84\x80y\xef\xfew\xaf\n\x03\x95h\xa1\xee\xda}D\xbf\x87\x10a\xf3\xc6\x92\xe7\xa3\xa3\x9dTR\tY'
    Na = 2
    authorities = list(range(Na))
    Nv_max = 10
    Nv = 6
    assert(Nv <= Nv_max)
    CA = {}
    BB = BBoard(CA, authorities)

    Nc = 5
    votes = [
        [2, 0, 0, 0, 0],
        [0, 1, 0, 0, 1],
        [0, 0, 1, 0, 0],
        [0, 1, 0, 0, 0],
        [0, 1, 0, 0, 0],
        [0, 1, 0, 0, 0],
    ]
    assert(len(votes) == Nv)
    PP = PublicParams(Na, Nc, Nv_max)
    emulate_voting(PP, BB, votes, public_seed)

    # reset_all()
    # Nc = 130
    # PP = PublicParams(Na, Nc, Nv_max)
    # votes = list(vote_from_candidate(PP.Nc, candidate) for candidate in [0, 1, 2, 128, 128, 129])
    # # creating invalid votes
    # votes[0][0] = 2
    # votes[1][0] = 1
    # emulate_voting(PP, BB, votes, public_seed)
