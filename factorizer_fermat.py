from math import isqrt
from factorizer_abstract import Factorizer, FactorList
from utility import extract_power_of_2

is_square = lambda n: n > 0 and isqrt(n) ** 2 == n

class FermatMethod(Factorizer):
    def factor(self):
        N = self.N
        factors = dict()
        if N % 2 == 0:
            k, N = extract_power_of_2(N)
            factors[2] = k
        x = isqrt(N)
        t = 2*x+1
        r = x*x-N
        while not is_square(r):
            r += t
            t += 2
        x = (t-1) // 2
        y = isqrt(r)
        factors[x-y] = 1
        factors[x+y] = 1
        return FactorList(factors)
