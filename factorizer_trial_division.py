from factorizer_abstract import Factorizer, FactorList
from utility import get_primes
from math import isqrt
from collections import defaultdict
from tqdm import tqdm

class TrialDivision(Factorizer):
    def factor(self, max_prime = None, verbose = False):
        process_iterator = tqdm if verbose else lambda x:x
        primes = get_primes(max_prime if max_prime else isqrt(self.N)+1)
        N = self.N
        factors = defaultdict(int)
        for p in process_iterator(primes):
            while N % p == 0:
                factors[p] += 1
                N //= p
            if N == 1:
                break
        if N != 1: factors[N] += 1
        return FactorList(dict(factors))


def test():
    factor1 = TrialDivision(600851475143).factor() # Project Euler Question
    assert factor1 == FactorList([71, 839, 1471, 6857])
    assert factor1.check_prime()
    from math import factorial
    factor2 = TrialDivision(factorial(25)**2).factor(max_prime = 100) # 25!^2
    assert factor2 == FactorList({
        2:44,
        3:20,
        5:12,
        7:6,
        11:4,
        13:2,
        17:2,
        19:2,
        23:2
    }) # 25!^2
    assert factor2.check_prime()
