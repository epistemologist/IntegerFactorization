from factorizer_abstract import Factorizer
from dataclasses import dataclass
from utility import Int, ModN, tonelli_shanks, get_primes, legendre, isqrt, count_primes
from typing import *
from itertools import count
from tqdm import tqdm

import numpy as np
import pdb

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
            if b != 0: return {int(-1*c/b)}
            else: return set()
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
            p: poly.mod_roots(p) for p in prime_base if poly.mod_roots(p) # Ignoring the 2 in prime_base
        }

class QuadraticSieve(Factorizer):

    def factor(self, use_mpqs = False, B = 100, chunk_size = 10**2):
        N = self.N
        prime_base = [2] + [p for p in get_primes(B) if legendre(N, p) == 1]
        """ Get values of x s.t. p(x) is smooth via a segmented sieve  """
        def _gen_candidates(factor_base: FactorBase, chunk_start: int, chunk_size: int, use_mpqs: bool):
            candidates = []
            # Sieve for smooth values
            sieve_chunk = np.zeros(shape = (chunk_size,))
            for p in factor_base.roots:
                # print(p, factor_base.roots[p])
                for root in factor_base.roots[p]:
                    sieve_chunk[(root - chunk_start) % p :: p] += np.log(p)
            # Check if vals are actually smooth and generate corresponding exponent vector
            for i in np.argpartition(sieve_chunk, -10)[-10:]:
                X = chunk_start + i
                val = factor_base.poly.eval(X)
                exponent_vector = np.zeros(shape = (len(prime_base), ), dtype=bool)
                tmp = val
                for i, p in enumerate(prime_base):
                    while tmp % p == 0:
                        tmp //= p
                        exponent_vector[i] += 1
                if tmp == 1:
                    if use_mpqs:
                        raise NotImplementedError()
                    else:
                        candidates.append((val, exponent_vector))
            return candidates
        def _poly_iterator():
            if not use_mpqs:
                poly = Polynomial(1, 0, -N)
                factor_base = FactorBase(poly, prime_base)
                for i in count(1):
                    yield (factor_base, isqrt(N) + i*chunk_size, chunk_size, False)
            else:
                raise NotImplementedError()
        smooth_vals = []
        for args in _poly_iterator():
            print(args, len(smooth_vals))
            smooth_vals.extend(_gen_candidates(*args))
            if len(smooth_vals) > len(prime_base):
                break
        return smooth_vals

from utility import get_semiprime

N = 1593004271

print(QuadraticSieve(N).factor())