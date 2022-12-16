#!/usr/bin/python3

from utility import get_primes, legendre, tonelli_shanks, isqrt
from factorizer_abstract import FactorList
from time import time
from Crypto.Util.number import getPrime
from collections import defaultdict
import numpy as np
from itertools import count

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


def benchmark_sieving(N, BLOCK_LENGTH):
    def _is_smooth(X):
        factorization = defaultdict(int)
        for p in primes:
            while X % p == 0:
                X //= p
                factorization[p] += 1
        return FactorList(dict(factorization)) if X == 1 else None
    B = isqrt(L(N)) + 1 # We can modify this later
    BLOCK_LENGTH = 10**8 // 5
    CANDIDATE_LENGTH = 1000
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
    def _process_chunk(chunk_start):
        sieve = np.zeros(BLOCK_LENGTH, dtype=np.float16)
        for p in primes:
            for root in square_roots[p]:
                sieve[ ((-chunk_start + root) % p)::p ] -= np.log(p)
            for root in square_root_p_squared[p]:
                sieve[((-chunk_start + root) % p*p) :: p*p] -= np.log(p)
        candidate_indices = np.argpartition(sieve, CANDIDATE_LENGTH)[:CANDIDATE_LENGTH]
        candidates_out = dict()
        for x in candidate_indices:
            x = int(x)
            candidate = (x+chunk_start)*(x+chunk_start) - N
            factorization = _is_smooth(candidate)
            if factorization is not None: candidates_out[x+chunk_start] = factorization.to_dict()
        return candidates_out
    start = time()
    chunks = 0
    candidates = dict()
    for chunk_start in count(isqrt(N)+1, BLOCK_LENGTH):
        candidates.update(_process_chunk(chunk_start))
        chunks += 1
        if len(candidates) >= 2*len(primes):
            break
    end = time()
    return end-start, chunks

# print(benchmark_sieving(getPrime(48) * getPrime(48), 10**7))


"""
f = open("out.log", "w")

for BITS in range(32, 72+4, 4):
    for BLOCK_SIZE in sum([[pow(10, i), 3*pow(10, i)] for i in range(6, 9)], []):
        for trial in range(5):
            out = f"{BITS}, {BLOCK_SIZE}, {trial}, {benchmark_sieving(getPrime(BITS) * getPrime(BITS), BLOCK_SIZE)}"
            print(out)
            f.write(out + "\n")
"""
N = 13241789903837055721
benchmark_sieving(N, 10**8)
