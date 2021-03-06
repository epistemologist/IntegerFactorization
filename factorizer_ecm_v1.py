from elliptic_curve import *
from factorizer_abstract import Factorizer, FactorList
from random import randint
from math import gcd, log
from utility import get_primes, InversionError
from tqdm import tqdm
from concurrent.futures import ThreadPoolExecutor, ProcessPoolExecutor, as_completed
from multiprocessing import cpu_count


B = 10000

def factor(N):
    while True:
        x = randint(0, N-1)
        y = randint(0, N-1)
        a = randint(0, N-1)
        b = (y**2 - x**3 - a*x) % N
        g = gcd(4*a**3 + 27*b**2, N)
        if g == N:
            continue
        elif g > 1:
            return g
        else:
            break
    curve = EllipticCurve( N, a, b )
    pt = EllipticCurveAffinePoint(curve, x, y)
    for i in range(2, B):
        try:
            pt = pt * i
        except InversionError as e:
            return gcd(e.n, e.p)

class ECM(Factorizer):
    def factor(self, B = 1000000, num_curves = 5000, verbose=True, threaded = False):
        N = self.N
        pbar = None
        # Construct power ladder
        prime_ladder = {p: int(log(B, p)) for p in get_primes(B)}
        if verbose:
            pbar = tqdm(total= sum(prime_ladder.values()) * num_curves)
        def _attempt_curve():
            while True:
                x = randint(0, N-1)
                y = randint(0, N-1)
                a = randint(0, N-1)
                b = (y**2 - x**3 - a*x) % N
                g = gcd(4*a**3 + 27*b**2, N)
                if g == N: continue
                elif g > 1:
                    return g
                else: break
            curve = EllipticCurve(N, a, b)
            pt = EllipticCurveAffinePoint(curve, x, y)
            #for i in range(2,B+1):
            for p,e in prime_ladder.items():
                for _ in range(e):
                    pbar.update(1)
                    try:
                        pt = pt * p
                    except InversionError as err:
                        return gcd(err.n, err.p)
        if threaded:
            with ThreadPoolExecutor(max_workers = cpu_count()) as executor:
                futures = [executor.submit(_attempt_curve) for i in range(num_curves)]
                for future in as_completed(futures):
                    tmp = future.result()
                    if tmp:
                        return FactorList([tmp, N//tmp])

        else:
            for i in tqdm(range(num_curves)):
                tmp = _attempt_curve()
                if tmp:
                    return FactorList([tmp, N//tmp])

N = 285806789255211805294914089243273
