import pdb
import numpy as np
from collections import defaultdict, OrderedDict
from tqdm import tqdm
from itertools import count
from typing import Optional
from galois import GF

from factorizer_abstract import Factorizer, FactorList
from utility import get_primes, legendre, tonelli_shanks, isqrt

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


def gaussian_elimination(M):
    m, n = M.shape
    h = k = 0
    while h < m and k < n:
        max( range(h, m+1), key = lambda i: M[i, k]  )

# Vanilla quadratic sieve
class QuadraticSieve(Factorizer):
    def factor(self):
        # Initialize
        N = self.N
        B = isqrt(L(N)) + 1 # We can modify this later
        primes = [p for p in get_primes(B) if p % 2 == 1 and legendre(N, p) == 1]
        square_roots = {}
        square_root_p_squared = {}
        for p in primes:
            square_root = tonelli_shanks(N, p)
            if square_root is None:
                raise ValueError(f"could not find square root mod {p}")
            else:
                square_roots[p] = (square_root, (-square_root) % p)
                square_root_p_squared[p] = sqrt_mod_p_squared(N, p)

        # Sieving

        BLOCK_LENGTH = 10**8
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
            for p in tqdm(primes):
                for root in square_roots[p]:
                    sieve[ ((-chunk_start + root) % p)::p ] -= np.log(p)
                for root in square_root_p_squared[p]:
                    sieve[((-chunk_start + root) % p*p) :: p*p] -= np.log(p)
            candidate_indices = np.argpartition(sieve, CANDIDATE_LENGTH)[:CANDIDATE_LENGTH]
            candidates_out = dict()
            for x in candidate_indices:
                candidate = (x+chunk_start)*(x+chunk_start) - N
                factorization = _is_smooth(candidate)
                if factorization is not None: candidates_out[candidate] = factorization
            return candidates_out

        candidates = dict()
        for chunk_start in count(isqrt(N)+1, BLOCK_LENGTH):
            candidates.update(_process_chunk(chunk_start))
            print(len(candidates))
            if len(candidates) >= len(primes) + 1:
                break

        # Linear algebra
        candidates = OrderedDict([(c, candidates[c]) for c in list(sorted(candidates))[:len(primes)+1]])
        exponent_vectors = []
        for c in candidates:
            factorization = candidates[c].to_dict()
            exponent_vectors.append([factorization[p] % 2 for p in primes ])
        return exponent_vectors


N = 8754660968220887821
# N = 1413409093
x = QuadraticSieve(N).factor()
