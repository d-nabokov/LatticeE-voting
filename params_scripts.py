def get_q(l, starting=2):
    Pr = Primes()
    q = starting
    while not (q - 1) % (4*l) == 2*l:
        q = Pr.next(q)
    return q


def get_zeta(l, q):
    root = None
    for val in range(2, q):
        if val**(2*l) % q == 1 and not val**(l) % q == 1:
            root = val
            break
    return root


def get_ring_params(l):
    q = get_q(l, Integer('111111111000000000000000000000', 2))
    q_bits = len(bin(q)[2:])
    zeta = get_zeta(l, q)
    print(f'q = {self.q}\nzeta = {self.zeta}\n')
    return q, q_bits, zeta
