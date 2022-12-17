from random import randint, seed
from tqdm import tqdm
from utility import Int, gcd
from factorizer_abstract import Factorizer, FactorList


class PollardRho(Factorizer):
    def factor(self, num_candidates=100, max_iter=1000000, verbose = True, seed_ = None):
        N = self.N
        if verbose: pbar = tqdm(total=num_candidates*max_iter)
        if seed_: seed(seed_)
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
                    if verbose: pbar.update(-pbar.n % max_iter)
                    break
                if g != 1:
                    return FactorList([g, N//g])
                if verbose: pbar.update(1)
