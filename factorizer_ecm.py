from collections import namedtuple
from math import log, gcd
from factorizer_abstract import Factorizer, FactorList
from utility import ModN, get_primes

# Some utility classes
Curve = namedtuple("Curve", ["A", "B", "C"])
ProjectivePoint = namedtuple("ProjectivePoint", ["X", "Z"])

# implementation from "Prime Numbers: A Computational Perspective"

# Equation 7.6
def add_h(curve: Curve, P1: ProjectivePoint, P2: ProjectivePoint, Pminus: ProjectivePoint) -> ProjectivePoint:
    (A,B,C), (X1, Z1), (X2, Z2), (Xminus, Zminus) = curve, P1, P2, Pminus
    if X1 == Z1 == 0:
        return (X2, Z2)
    if X2 == Z2 == 0:
        return (X1, Z1)
    tmp = X1*X2 - A*Z1*Z2
    Xplus = Zminus*(tmp*tmp - 4*B*Z1*Z2*(X1*Z2 + X2*Z1 + C*Z1*Z2))
    tmp2 = X1*Z2 - X2*Z1
    Zplus = Xminus*tmp2*tmp2
    return ProjectivePoint(Xplus, Zplus)

# Equation 7.7

def double_h(curve: Curve, P: ProjectivePoint):
    (A,B,C), (X1, Z1) = curve, P
    tmp = X1*X1 - A*Z1*Z1
    Xplus = tmp*tmp - 4*B*Z1*Z1*Z1*(2*X1 + C*Z1)
    Zplus = 4*Z1*(X1*X1*X1 + C*X1*X1*Z1 + A*X1*Z1*Z1 + B*Z1*Z1*Z1)
    return ProjectivePoint(Xplus, Zplus)

# Algorithm 7.2.7
def point_multiplication(curve: Curve, P: ProjectivePoint, n: int) -> ProjectivePoint:
    (A,B,C), (X,Z) = curve, P
    # if n == 0: return [0, 0]
    if n == 0: return ProjectivePoint(0,0)
    if n == 1: return P
    if n == 2: return double_h(curve, P)
    U = P
    T = double_h(curve, P)
    bit_length = int(n).bit_length()
    for j in range(bit_length - 2, 0, -1):
        if (n >> j) & 1:
            U = add_h(curve, T, U, P)
            T = double_h(curve, T)
        else:
            T = add_h(curve, U, T, P)
            U = double_h(curve, U)
    if n & 1: return add_h(curve, U,T,P)
    return double_h(curve, U)

class ECM(Factorizer):
    def factor(self, B1 = 10000, B2 = 50000, D = 100, num_curves = 100):
        assert B1 % 2 == 0 and B2 % 2 == 0
        N = self.N
        while N % 2 == 0:
            N //= 2
        while N % 3 == 0:
            N //= 6
        power_ladder = { p: int(log(B1, p)) for p in get_primes(B1)} # Used in stage 1
        primes = get_primes(B2 + 2*D) # Used in stage 2
        prime_map = {r: [q for q in primes if r+2 <= q <= r+2*D] for r in range(B1-1, B2, 2*D) }
        # Implementation of Algorithm 7.4.4
        def _process_curve(seed): # Here, seed is sigma
            u = ModN(seed*seed - 5, N)
            v = ModN(4*seed, N)
            C = ((v-u)**3 * (3*u+v))/(4*u**3*v) - 2 # Underlying curve is y^2 = x^3 + Cx^2 + x
            # i.e. (A,B,C) = (1,0,C)
            curve = Curve(1,0,C)
            Q = ProjectivePoint(u**3, v**3)
            # Stage 1
            for p, e in power_ladder.items():
                Q = point_multiplication(curve, Q, p**e)
            g = gcd(int(Q.Z), N)
            if 1 < g < N: return g
            # Stage 2
            S1 = double_h(Q)
            S2 = double_h(S1)
            S = [None, S1, S2]
            beta = [None]
            for d in range(1, D+1):
                if d > 2:
                    S.append(add_h(C, S[d-1], S[1], S[d-2]))
                curr_S = S[d]
                beta.append( curr_S.X * curr_S*Z )
            g = ModN(1, N)
            B = B1 - 1
            T = point_multiplication(curve, Q, B - 2*D)
            R = point_multiplication(curve, B, Q)
            for r in prime_map:
                alpha = R.X * R.Z
                for q in prime_map[r]:
                    delta = (q-r) // 2
                    S_delta = S[delta]
                    # g = g((X(R) - X(S_d))(Z(R) + Z(S_d))-a+b_d)
                    g *= (R.X - S_delta.X) * (R.Z + S_delta.Z) - alpha + beta[delta]
                R, T =add_h(R, S[D], T),R
            g = gcd(int(g),n)
            if 1 < g < n: return g
            else: return None
        for i in range(6, num_curves+6):
            temp = _process_curve(i)
            if temp: return temp
