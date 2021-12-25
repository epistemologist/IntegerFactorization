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

def is_perfect_power(N: int, n: int) -> bool:
     # Returns whether or not there exists an x s.t. N = x^n
     # Extract kth root of N using binary search
     def _integer_root(n,k, max_iter = 10000):
          hi = 1
          while pow(hi, k) < n: hi *= 2
          lo = hi // 2
          num_iter = 0
          while num_iter < max_iter:
              mid = (lo + hi) // 2
              tmp = pow(mid, k)
              if tmp == n:
                  return (mid, mid)
              if hi-lo == 1:
                  break
              num_iter += 1
              if tmp < n:
                  lo = mid
              else:
                  hi = mid
          return (lo, hi)
     lo, hi = _integer_root(N, n)
     return pow(lo, n) == N or pow(hi, n) == N

# Simple utility class to represent integers mod prime
class ModN:
    def __init__(self, n: int, p: int):
        self.n = n
        self.p = p
    def _shares_modulus(self, other):
        return isinstance(other, ModN) and self.p == other.p
    def __eq__(self, other):
        if type(other) == int: other = ModN(other, self.n)
        return self._shares_modulus(other) and self.n % self.p == other.n % other.p
    def __pow__(self, other):
        if type(other) != int: raise ValueError("cannot raise number to non-integer exponent")
        return ModN(pow(self.n, other, self.p), self.p)
    def __repr__(self):
        return str(self.n)

def _create_internal_fn(self, other, func):
    if type(other) == int: other = ModN(other, self.p)
    if not self._shares_modulus(other): raise ValueError("Cannot do operation!")
    return func(self, other)

ModN.__add__ = lambda x,y: _create_internal_fn(x,y,lambda s,o: ModN((s.n + o.n) % s.p,s.p))
ModN.__radd__ = lambda x,y: _create_internal_fn(x,y,lambda s,o: ModN((o.n + s.n) % s.p,s.p))
ModN.__sub__ = lambda x,y: _create_internal_fn(x,y,lambda s,o: ModN((s.n - o.n) % s.p,s.p))
ModN.__rsub__ = lambda x,y: _create_internal_fn(x,y,lambda s,o: ModN((o.n - s.n) % s.p,s.p))
ModN.__mul__ = lambda x,y: _create_internal_fn(x,y,lambda s,o: ModN((s.n * o.n) % s.p,s.p))
ModN.__rmul__ = lambda x,y: _create_internal_fn(x,y,lambda s,o: ModN((o.n * s.n) % s.p,s.p))
ModN.__truediv__ = lambda x,y: _create_internal_fn(x,y,lambda s,o: ModN((s.n * pow(o.n, -1, s.p)) % s.p,s.p))
ModN.__rtruediv__ = lambda x,y: _create_internal_fn(x,y,lambda s,o: ModN((o.n * pow(s.n, -1, s.p) % s.p,s.p)))
