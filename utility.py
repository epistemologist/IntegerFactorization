from math import log
from typing import List, Tuple
from random import randint
from tqdm import tqdm
import prime_sieve as sieve

def get_primes(N: int) -> List[int]:
    # Returns all primes p such that 1 < p <= N
    if N < 10**6:
        return sieve.get_small_primes(N)
    else:
        N_upper = (int(N**0.5)+1)**2
        s = sieve.FastPrimeSieve(N_upper, threading=False, progress_bar=False)
        s.get_primes()
        return [i for i in s.primes if i <= N]

def extract_power_of_2(N: int) -> Tuple[int, int]:
    # Given N, this function returns integers k,d such that 2^k * d == N
    k = 0
    d = N
    while d % 2 == 0:
        d //= 2
        k += 1
    return (k,d)

def naive_is_prime(N: int) -> bool:
    # Naive primality test by trial division
    if N < 2: return False
    if N == 2: return True
    for d in range(2, int(N**0.5) + 1):
        if N % d == 0: return False
    return True

def miller_rabin(N: int, k: int) -> bool:
    # Miller Rabin test
    if k == 0: # deterministic
        witnesses = range(2, min(N-1,int(2*log(N)**2)))
    else: # probabilistic
        witnesses = map(lambda _: randint(2, N-2), range(k))
    r,d = extract_power_of_2(N-1)
    for a in witnesses:
        x = pow(a,d,N)
        if x == 1 or x == N-1: continue
        for _ in range(r-1):
            x = pow(x,2,N)
            if x == N-1:
                break
        else:
            return False
    return True

def is_prime(N: int) -> bool:
    if N < 10**5:
        return naive_is_prime(N)
    else:
        return miller_rabin(N, 50)
