from factorizer_abstract import Factorizer
from dataclasses import dataclass
from utility import Int, ModN, tonelli_shanks, get_primes, legendre, isqrt
from typing import *
from itertools import count
from tqdm import tqdm

import numpy as np

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
    def __init__(self, poly: Polynomial, B: int):
        self.poly = poly
        self.roots = {
            p: poly.mod_roots(p) for p in get_primes(B) if poly.mod_roots(p)
        }

class QuadraticSieve(Factorizer):
    """ Get potentially smooth values of a polynomial via a segmented sieve  """
    def _get_smooth_vals(self, factor_base: FactorBase, chunk_start: int, chunk_size: int):
        sieve_chunk = np.zeros(shape = (chunk_size,))
        for p in tqdm(factor_base.roots):
            for root in factor_base.roots[p]:
                sieve_chunk[(root - chunk_start) % p :: p] += np.log(p)
        smooth_vals = [factor_base.poly.eval(chunk_start + i) for i in np.argpartition(sieve_chunk, -10)[-10:]]
        return smooth_vals

    def factor(self, use_mpqs = False, B = 10**6, chunk_size = 10**5):
        N = self.N
        def _poly_iterator():
            if not use_mpqs:
                poly = Polynomial(1, 0, -N)
                factor_base = FactorBase(poly, B)
                for i in count(1):
                    yield (factor_base, isqrt(N) + i*chunk_size, chunk_size)
            else:
                ...
        smooth_vals = []
        for args in _poly_iterator():
            smooth_vals.extend(self._get_smooth_vals(*args))
            if len(smooth_vals) > 100:
                break
        return smooth_vals
