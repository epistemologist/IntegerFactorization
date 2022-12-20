from utility import get_semiprime
from factorizer_abstract import Factorizer
from dataclasses import dataclass
from utility import Int, ModN, tonelli_shanks, get_primes, legendre, isqrt, gcd
from typing import *
from itertools import count
from tqdm import tqdm
from math import prod, sqrt
from multiprocessing import Pool


import numpy as np

import galois
import pdb

VERBOSE_LEVEL = None

""" Used to calculate best B corresponding to n """
def L(n: Int) -> Int:
    ln_n = n.bit_length()
    ln_ln_n = ln_n.bit_length()
    return Int(1) << Int(sqrt(ln_n * ln_ln_n))


@dataclass
class Polynomial:
    a: Int
    b: Int
    c: Int

    def eval(self, x: Int) -> Int:
        return self.a*x*x + self.b*x + self.c

    def mod_roots(self, p: Int) -> Set[Int]:
        a, b, c = [ModN(i, p) for i in (self.a, self.b, self.c)]
        if a == 0:
            if b != 0:
                return {int(-1*c/b)}
            else:
                return set()
        else:
            D = b*b - 4*a*c
            try:
                d = D.sqrt()
            except:
                return None
            r1 = (-1*b + d) / (2*a)
            r2 = (-1*b - d) / (2*a)
            return {int(r1), int(r2)}


""" Factor base utility class used to hold information for sieving  """


class FactorBase:
    def __init__(self, poly: Polynomial, prime_base: List[int]):
        self.poly = poly
        self.roots = {
            p: poly.mod_roots(p) for p in prime_base
        }


"""Given matrix M over F_2, give a list of its rows s.t. sum of rows equals 0"""


def get_null_space(M: np.array) -> List[List[int]]:
    # Default to using Galois for right now
    # TODO: Implement faster algorithm
    return [[n for n, i in enumerate(null_space_vector) if i] for null_space_vector in galois.GF2(M).T.null_space()]


"""Returns corresponding exponent vector for N if N is smooth w.r.t given prime base"""


def _get_exponent_vector(N: Int, prime_base: List[int]) -> Optional[np.array]:
    tmp = N
    exponent_vector = np.zeros(shape=(len(prime_base),), dtype=int)
    for i, p in enumerate(prime_base):
        while tmp % p == 0:
            tmp //= p
            exponent_vector[i] += 1
    return None if tmp != 1 else exponent_vector


""" Get values of x s.t. p(x) is smooth via a segmented sieve  """


def _gen_candidates(factor_base: FactorBase, chunk_start: int, chunk_size: int, use_mpqs: bool):
    """
    Inputs:
        factor_base: FactorBase object that contains primes as well as polynomial to sieve
        chunk_start: value of X to start sieving at
        chunk_size: size of chunk to sieve
        use_mpqs: whether or not we are using MPQS (defaults to using naive single polynomial for sieving)
    Output:
        candidates - list of [(X, p(X), exponent vector of p(X)) ] such that p(X) is smooth w.r.t given factor base
    """
    global VERBOSE_LEVEL
    candidates = []
    # Sieve for smooth values
    sieve_chunk = np.zeros(shape=(chunk_size,))
    for p in factor_base.roots:
        # print(p, factor_base.roots[p])
        if p == 2:
            continue
        if factor_base.roots[p]:
            for root in factor_base.roots[p]:
                sieve_chunk[(root - chunk_start) % p:: p] += np.log(p)
    # Check if vals are actually smooth and generate corresponding exponent vector
    for i in np.argpartition(sieve_chunk, -10)[-10:]:
        X = chunk_start + i
        val = factor_base.poly.eval(X)
        exponent_vector = _get_exponent_vector(val, list(factor_base.roots))
        if exponent_vector is not None:
            if use_mpqs:
                raise NotImplementedError()
            else:
                candidates.append((X, val, exponent_vector))
    return candidates


def _sieve_worker(args):
    global VERBOSE_LEVEL
    if VERBOSE_LEVEL >= 3: print(args)
    return _gen_candidates(*args)


class QuadraticSieve(Factorizer):
    def factor(self, use_mpqs=False, B=None, chunk_size=10**6, num_processes=2, verbosity=2):
        """
        Inputs:
            use_mpqs: if True, will use MPQS instead of single polynomial (not implemented)
            B: smoothness bound (if left uninitialized, will be automatically calculated based on size of N)
            chunk_size: size of chunk to be sieved by each worker
            num_proceses: number of cores to dedicate to sieving
            verbosity:
             - 0: supress all output
             - 1: print different stages
             - 2: have progress bars for sieving and linear algebra
             - 3: print everything
        Output:
            factorization
        """
        N = self.N
        if B is None:
            B = isqrt(L(N))
        prime_base = [2] + \
            [p for p in get_primes(B) if p > 2 and legendre(N, p) == 1]
        global VERBOSE_LEVEL
        VERBOSE_LEVEL = verbosity
        """ Infinite iterator to generate args for sieve workers """
        def _poly_iterator():
            if not use_mpqs:
                poly = Polynomial(1, 0, -N)
                factor_base = FactorBase(poly, prime_base)
                for i in count(1):
                    yield (factor_base, isqrt(N) + i*chunk_size, chunk_size, False)
            else:
                raise NotImplementedError()
        # Generate smooth candidates
        candidates = []
        candidates_pbar = tqdm(total = len(prime_base)) if VERBOSE_LEVEL == 2 else None
        with Pool(processes=num_processes) as p:
            for result in p.imap_unordered(_sieve_worker, _poly_iterator()):
                candidates.extend(result)
                if candidates_pbar is not None: candidates_pbar.update(len(result))
                if VERBOSE_LEVEL >= 3: print(f"[*] Collected {len(candidates)} out of {len(prime_base)} candidates")
                if len(candidates) >= len(prime_base):
                    p.terminate()
                    break
        if VERBOSE_LEVEL >= 1:
            print("[!] Completed sieving")
        # Linear algebra
        M = np.vstack(i[2] for i in candidates).astype(np.uint32) % 2
        null_space_vectors = get_null_space(M)
        if VERBOSE_LEVEL >= 1:
            print("[!] Completed linear algebra")
        for v in null_space_vectors:
            new_candidates = [candidates[i] for i in v]
            # Factor
            x = prod([i[0] for i in new_candidates])
            y = prod([Int(p**e) for p, e in zip(prime_base,
                                                sum([i[2] for i in new_candidates]) // 2)])
            d = gcd(x-y, N)
            if d != 1 and d != N:
                return d
        return None

N = 16403383910217002297

out = QuadraticSieve(N).factor()
