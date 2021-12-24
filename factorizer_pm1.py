from math import log, gcd, prod
from bisect import bisect
from factorizer_abstract import Factorizer, FactorList
from utility import get_primes

lcm = lambda x,y: (x*y)//gcd(x,y)
LCM = lambda arr: reduce(lcm, arr)
from functools import reduce

class PM1(Factorizer):
    def factor(self, B1 = 10000, B2=25000):
        assert B1 < B2
        N = self.N
        primes = get_primes(B2)
        divider = bisect(primes, B1)
        small_primes, large_primes = primes[:divider], primes[divider:]
        # Stage 1
        for a in [2,3,5,7]:
            # Calculate a_Q = a^Q % N, Q = \prod_{primes q <= B} q^(int(log(B,q)))
            a_Q = a
            for q in small_primes:
                for i in range(int(log(B1,q))):
                    a_Q = pow(a_Q, q, N)
            g = gcd(a_Q-1, N)
            print(a, g, a_Q, pow(a, LCM(range(1, B1+1)),N))
            if g in [1, N]:
                # Stage 2
                prime_differences = [j-i for i,j in zip(large_primes, large_primes[1:])]
                # Calculate {a^(Q*q_i) for q_i in large_primes} by using differences
                # If we have a^(Q*q_k), we can calculate a^(Q*q_{k+1})
                # a^(Q*q_{k+1}) = a^(Q*[q_k + (q_{k+1} - q_k)])
                # = a^(Q*q_k) * a^(Q * [q_{k+1} - q_k] )
                # we can precompute the second factor
                cache = dict() # cache[i] = (a^(Q*i))
                tmp = 1
                for i in range(0, max(prime_differences)+2, 2):
                    cache[i] = tmp
                    tmp = (tmp * a_Q * a_Q) % N
                pow_a_large_prime = pow(a_Q, large_primes[0], N)
                for q_1, q_2 in zip(large_primes, large_primes[1:]):
                    g = gcd((pow_a_large_prime - 1) % N, N)
                    if g not in [1, N]:
                        return FactorList([g, N//g])
                    pow_a_large_prime *= cache[q_2-q_1]
                    pow_a_large_prime %= N
                    # print(pow_a_large_prime, pow(a_Q, q_2, N))
            else:
                return FactorList([g, N//g])

