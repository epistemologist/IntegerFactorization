import pdb
import numpy as np
import bottleneck as bn
from collections import defaultdict, OrderedDict
from tqdm import tqdm
from itertools import count
from typing import Optional, List
from math import prod, gcd
from random import sample
from numba import njit
from time import time
from Crypto.Util.number import getPrime


from factorizer_abstract import Factorizer, FactorList
from utility import get_primes, legendre, tonelli_shanks, isqrt, Timer

def L(n):
    ln_n = n.bit_length()
    ln_ln_n = ln_n.bit_length()
    return pow(2, isqrt(ln_n * ln_ln_n))

def sqrt_mod_p_squared(n, p):
    # Based on https://mathoverflow.net/a/223806
    x = tonelli_shanks(n, p)
    if x is None:
        return None
    return [
        (pow(root, p, p*p) * pow(n, (p*p-p+1) // 2, p*p)) % (p*p)
        for root in [x, -x%p]
    ]


# Vanilla quadratic sieve
class QuadraticSieve(Factorizer):
    def factor(self, VERBOSE = 3, B = None):
        """
        Verbosity parameter:
        VERBOSE = 0: print nothing
        VERBOSE = 1: print breakdown of main stages
        VERBOSE = 2: print parameters
        VERBOSE = 3: have progress bars
        """
        progress_bar = lambda args: args if VERBOSE < 3 else tqdm(args)
        # Initialize
        N = self.N
        if VERBOSE > 1: print(N)
        if VERBOSE > 0: print(f"Starting and initializing parameters... "); timer = Timer()
        if B is None: B = 2 * (isqrt(L(N)) + 1) # We can modify this later
        primes = [p for p in get_primes(B) if p % 2 == 1 and legendre(N, p) == 1]
        if VERBOSE > 1: print(f"primes: {min(primes), max(primes)}")
        square_roots = {}
        square_root_p_squared = {}
        for p in primes:
            square_root = tonelli_shanks(N, p)
            if square_root is None:
                raise ValueError(f"could not find square root mod {p}")
            else:
                square_roots[p] = (square_root, (-square_root) % p)
                square_root_p_squared[p] = sqrt_mod_p_squared(N, p)
        if VERBOSE > 0: print(f"Finished finding square roots: {timer()}")
        # Sieving

        BLOCK_LENGTH = 10**8 // 5
        CANDIDATE_LENGTH = 1000

        def _is_smooth(X) -> Optional[FactorList]:
            factorization = defaultdict(int)
            for p in primes:
                while X % p == 0:
                    X //= p
                    factorization[p] += 1
            return FactorList(dict(factorization)) if X == 1 else None

        def _process_chunk(chunk_start):
            sieve = np.zeros(BLOCK_LENGTH, dtype=np.float16)
            for p in progress_bar(primes):
                for root in square_roots[p]:
                    sieve[ ((-chunk_start + root) % p)::p ] -= np.log(p)
                for root in square_root_p_squared[p]:
                    sieve[((-chunk_start + root) % p*p) :: p*p] -= np.log(p)
            candidate_indices = bn.argpartition(sieve, CANDIDATE_LENGTH)[:CANDIDATE_LENGTH]
            candidates_out = dict()
            for x in candidate_indices:
                x = int(x)
                candidate = (x+chunk_start)*(x+chunk_start) - N
                factorization = _is_smooth(candidate)
                if factorization is not None: candidates_out[x+chunk_start] = factorization.to_dict()
            return candidates_out

        candidates = dict() # dict of {x: x^2 - N} s.t. x^2-N is smooth
        for chunk_start in count(isqrt(N)+1, BLOCK_LENGTH):
            candidates.update(_process_chunk(chunk_start))
            if VERBOSE > 0: print(f"Processed chunk with chunk_start {chunk_start}: {timer()}")
            if VERBOSE > 1: print(len(candidates), 2*len(primes))
            return
            if len(candidates) >= 2*len(primes):
                break
        if VERBOSE > 0: print(f"Finished sieving: {timer()}")
        # While we haven't found a factor...
        while True:
            # Linear algebra
            candidates = OrderedDict([(c, candidates[c]) for c in sample(list(candidates), len(primes)+1)])
            exponent_vectors = []
            for c in candidates:
                factorization = candidates[c]
                exponent_vectors.append([factorization[p] % 2 for p in primes ])
            # Inspired by https://github.com/mikolajsawicki/quadratic-sieve/blob/main/quadratic_sieve/fast_gauss.py
            # Also referred to https://www.cs.umd.edu/~gasarch/TOPICS/factoring/fastgauss.pdf
            M = np.vstack(exponent_vectors)
            if VERBOSE > 1: print(f"M: {M}, {M.shape}")
            @njit()
            def _to_rref(M):
                n,m = M.shape
                marked = []
                for j in range(m):
                    if j % max(100, m // 100) == 0 and VERBOSE > 2: print(j, m)
                    i = 0
                    while i < n and M[i,j] != 1:
                        i += 1
                    if i == n:
                        continue
                    else:
                        marked.append(i)
                        for k in range(m):
                            if k != j:
                                if M[i, k] == 1:
                                    M[:, k] ^= M[:, j]
                return marked
            marked = _to_rref(M)
            if VERBOSE > 1: print(f"marked: {marked}")
            dependent_index = [i for i in range(len(M)) if i not in marked][0]
            if VERBOSE > 1: print(f"dependent index: {dependent_index}")
            basis = [dependent_index]
            ones_cols = np.where(M[dependent_index] == 1)[0]
            if VERBOSE > 1: print(f"ones_cols: {ones_cols}")
            basis.extend([ np.where(M[:, j] == 1)[0][0] for j in ones_cols ])
            if VERBOSE > 1: print(f"basis: {basis}")
            # Factorization step
            x_s = []
            y_s = []
            C = list(candidates.keys())
            for b in basis:
                x_s.append(C[b])
                y_s.append(candidates[C[b]])
            def _sqrt_prod( nums: List[FactorList] ) -> FactorList:
                out = defaultdict(int)
                for n in nums:
                    for p in n:
                        out[p] += n[p]
                if any([i % 2 != 0 for i in out.values()]): raise ValueError("not perfect square!")
                return FactorList({ p: (out[p]) // 2 for p in out })
            x = prod(x_s)
            y = _sqrt_prod(y_s).prod()
            d = gcd(x-y, N)
            if d != 1 and d != N:
                if VERBOSE > 0: print(f"Finished linear algebra: {timer()}")
                return d, N//d
            else:
                if VERBOSE > 1: print("Found trivial factorization, trying again...")

"""
from Crypto.Util.number import getPrime
import cProfile
k = 64
N = getPrime(k) * getPrime(k)

print(N)
"""

N = 162150089453462133245257817849249542687


for i in range(100):
    x = QuadraticSieve(N).factor(VERBOSE=1)



