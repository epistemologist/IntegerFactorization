from math import log, sqrt
from typing import List, Tuple
from random import randint
import prime_sieve as sieve

try:
    import gmpy2
    GMPY_IMPORT = True
    Int = lambda n: gmpy2.mpz(int(n))
except:
    GMPY_IMPORT = False
    Int = int

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

def legendre(N: int, p: int) -> int:
    # Assuming p is prime
    if N % p == 0: return 0
    # Using Euler's criterion
    N %= p
    N_over_p = pow(N, (p-1) // 2, p)
    return 1 if N_over_p == 1 else -1

def tonelli_shanks(N: int, p: int) -> int:
    # Assuming p is prime
    # returns x s.t. x^2 = N (mod p)
    # implementation of Algorithm 2.3.8 from Prime Numbers: A Computational Perspective
    N %= p
    if legendre(N, p) != 1:
        return None
    if p % 8 == 3 or p % 8 == 7:
        return pow(N, (p+1)//4, p)
    elif p % 8 == 5:
        x = pow(N, (p+3) // 8, p)
        c = (x*x) % p
        if (c % p != N % p):
            x *= pow(2, (p-1)//4, p)
            x %= p
        return x
    else:
        d = 2
        while legendre(d, p) != -1:
            d = randint(2, p-1)
        s, t = extract_power_of_2(p-1)
        A = pow(N,t,p)
        D = pow(d,t,p)
        m = 0
        for i in range(s):
            if pow(A * pow(D,m,p), pow(2, s-1-i), p) == p-1:
                m += pow(2, i)
        x = pow(N, (t+1)//2, p) * pow(D, m // 2, p)
        x %= p
        return x


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

def isqrt(N):
    if not GMPY_IMPORT:
        from math import isqrt
        return isqrt(N)
    else:
        from gmpy2 import isqrt
        return isqrt(N)


# Simple utility class to represent integers mod prime

# Custom exception to raise in case of failed inversion (for ECM factoring)

class InversionError(Exception):
    def __init__(self, n, p):
        self.n = n
        self.p = p
        super().__init__(f"given n ({n}) is not invertible given base {p}")

def inv(n, p):
    if not GMPY_IMPORT:
        try:
            return pow(n, -1, p)
        except:
            raise InversionError(n, p)
    else:
        import gmpy2
        try:
            return gmpy2.invert(n, p)
        except:
            raise InversionError(n, p)

def _create_internal_fn(func):
    def wrapped(self, other):
        if type(other) != ModN: other = ModN(other, self.p)
        if not self._shares_modulus(other):
            print(type(self), type(other))
            raise ValueError("Cannot do operation!")
        return func(self, other)
    return wrapped


class ModN:
    def __init__(self, n: int, p: int):
        self.n = Int(n) % Int(p)
        self.p = Int(p)
    def _shares_modulus(self, other):
        return isinstance(other, ModN) and self.p == other.p
    def __eq__(self, other):
        if type(other) == int: other = ModN(other, self.p)
        return self._shares_modulus(other) and self.n % self.p == other.n % other.p
    def __pow__(self, other):
        if type(other) != int: raise ValueError("cannot raise number to non-integer exponent")
        return ModN(pow(self.n, other, self.p), self.p)
    def __repr__(self):
        return str(self.n)
    def __hash__(self):
        return hash(self.n)
    def __int__(self):
        return int(self.n)

ModN.__add__ = _create_internal_fn(lambda s,o: ModN((s.n + o.n) % s.p,s.p))
ModN.__radd__ = _create_internal_fn(lambda s,o: ModN((o.n + s.n) % s.p,s.p))
ModN.__sub__ = _create_internal_fn(lambda s,o: ModN((s.n - o.n) % s.p,s.p))
ModN.__rsub__ = _create_internal_fn(lambda s,o: ModN((o.n - s.n) % s.p,s.p))
ModN.__mul__ = _create_internal_fn(lambda s,o: ModN((s.n * o.n) % s.p,s.p))
ModN.__rmul__ = _create_internal_fn(lambda s,o: ModN((o.n * s.n) % s.p,s.p))
ModN.__truediv__ = _create_internal_fn(lambda s,o: ModN((s.n * inv(o.n, s.p)) % s.p,s.p))
ModN.__rtruediv__ = _create_internal_fn(lambda s,o: ModN((o.n * inv(s.n, s.p) % s.p,s.p)))
