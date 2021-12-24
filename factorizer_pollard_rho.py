from random import randint
from math import gcd
from tqdm import tqdm
from gmpy2 import mpz
from factorizer_abstract import Factorizer, FactorList

class PollardRho(Factorizer):
    def factor(self,num_candidates=100 ,max_iter = 1000000, verbose = True):
        N = self.N
        if verbose:
            pbar = tqdm(total=num_candidates * max_iter)
            total_iter = 0
        for i in range(num_candidates):
            a = randint(1, N-3)
            s = randint(0, N-1)
            U = V = s
            F = lambda x: (x*x + a) % N
            for j in range(max_iter):
                U = F(U)
                V = F(F(V))
                g = gcd(abs(U-V), N)
                if g == N:
                    if verbose:
                        pbar.update(-total_iter % max_iter)
                        total_iter += -total_iter % max_iter
                    break
                if g != 1:
                    return FactorList([g, N//g])
                if verbose:
                    total_iter += 1
                    pbar.update(1)
