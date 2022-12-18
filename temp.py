from factorizer_qs_v2 import *
from utility import get_semiprime, get_primes, isqrt
from random import randint
from factorizer_pollard_rho import PollardRho

N = get_semiprime(48)
B = 10**6


p1 = Polynomial(1, 0, -N)
factor_base1 = {p for p in get_primes(B) if p1.mod_roots(p)}

# Polynomial f(x) = a x^2 + 2 b x + c s.t. b^2 - a c = n

b = (isqrt(N) + 100)
a, c = list(PollardRho(b**2 - N).factor().factor_list)

p2 = Polynomial(a, 2*b, c)
factor_base2 = {p for p in get_primes(B) if p2.mod_roots(p)}

print(factor_base1 == factor_base2)