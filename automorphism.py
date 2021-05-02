

def phi(PP, x, power):
    if power < 0:
        power += PP.k
    phi_index = 2 * PP.l // PP.k + 1
    res = 0
    for j, coef in enumerate(x):
        p = j * phi_index**power
        sign = (-1)**((p // PP.d) % 2)
        newpower = p % PP.d
        res += coef * sign * PP.X**(newpower)
    return res
    

if __name__ == '__main__':
    from params import PublicParams
    PP = PublicParams(2, 127, 10)
    X = PP.X
    
    assert(X**3 == phi(PP, phi(PP, X**3, 2), -2))
    assert(X**3 == phi(PP, phi(PP, phi(PP, X**3, 1), -3), 2))
    from random_polynomials import random_poly
    from utils import randombytes
    seed = randombytes(PP.seedlen)
    poly, _ = random_poly(PP, seed, 0, PP.d)
    assert(poly == phi(PP, phi(PP, poly, 1), -1))
    assert(poly == phi(PP, phi(PP, phi(PP, poly, 1), -3), 2))

