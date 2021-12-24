from factorizer_abstract import Factorizer, FactorList
from math import isqrt, gcd
# Helper functions
is_square = lambda N: isqrt(N) ** 2 == N
ceil_sqrt = lambda N: isqrt(N) if is_square(N) else isqrt(N)+1
class HartMethod(Factorizer):
    def factor(self):
        N=self.N;s=next(filter(lambda s:is_square(pow(s,2,N)), map(lambda i: ceil_sqrt(N*i), range(1,N+1))));t=isqrt(pow(s,2,N));f1=gcd(s-t,N)
        """
        N = self.N;
        s=next(filter(
            lambda s: is_square(pow(s,2,N)),
            map(
                lambda i: ceil_sqrt(N*i),
                range(1, N+1)
            )
        ));
        t=isqrt(pow(s,2,N));
        f1=gcd(s-t,N);
        """
        return FactorList([f1, N//f1])
