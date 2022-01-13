from collections import namedtuple
from math import log, gcd
from factorizer_abstract import Factorizer, FactorList
from utility import ModN, get_primes

# Some utility classes
Curve = namedtuple("Curve", ["A", "B", "C"])
ProjectivePoint = namedtuple("ProjectivePoint", ["X", "Z"])

# implementation from "Prime Numbers: A Computational Perspective"

# Equation 7.6
def add_h(A, B, C, X1, Z1, X2, Z2, Xminus, Zminus):
    if X1 == Z1 == 0:
        return (X2, Z2)
    if X2 == Z2 == 0:
        return (X1, Z1)
    tmp = X1*X2 - A*Z1*Z2
    Xplus = Zminus*(tmp*tmp - 4*B*Z1*Z2*(X1*Z2 + X2*Z1 + C*Z1*Z2))
    tmp2 = X1*Z2 - X2*Z1
    Zplus = Xminus*tmp2*tmp2
    return (Xplus, Zplus)

def add_h2(curve: Curve, P1: ProjectivePoint, P2: ProjectivePoint, Pminus: ProjectivePoint) -> ProjectivePoint:
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
def double_h(A, B, C, X1, Z1):
    tmp = X1*X1 - A*Z1*Z1
    Xplus = tmp*tmp - 4*B*Z1*Z1*Z1*(2*X1 + C*Z1)
    Zplus = 4*Z1*(X1*X1*X1 + C*X1*X1*Z1 + A*X1*Z1*Z1 + B*Z1*Z1*Z1)
    return (Xplus, Zplus)

def double_h(curve: Curve, P: ProjectivePoint):
    (A,B,C), (X1, Z1) = curve, P
    tmp = X1*X1 - A*Z1*Z1
    Xplus = tmp*tmp - 4*B*Z1*Z1*Z1*(2*X1 + C*Z1)
    Zplus = 4*Z1*(X1*X1*X1 + C*X1*X1*Z1 + A*X1*Z1*Z1 + B*Z1*Z1*Z1)
    return ProjectivePoint(Xplus, Zplus)

# Algorithm 7.2.7
def point_multiplication(A, B, C, X, Z, n):
    if n == 0: return [0, 0]
    if n == 1: return [X, Z]
    if n == 2: return double_h(A, B, C, X, Z)
    U, V = X, Z
    T, W = double_h(A, B, C, X, Z)
    bit_length = int(n).bit_length()
    for j in range(bit_length - 2, 0, -1):
        if (n >> j) & 1:
            U, V = add_h(A, B, C, T, W, U, V, X, Z)
            T, W = double_h(A, B, C, T, W)
        else:
            T, W = add_h(A, B, C, U, V, T, W, X, Z)
            U, V = double_h(A, B, C, U, V)
    if n & 1: return add_h(A, B, C, U, V, T, W, X, Z)
    return double_h(A, B, C, U, V)

def point_multiplication2(curve: Curve, P: ProjectivePoint, n: int):
	(A,B,C), (X,Z) = curve, P
    if n == 0: return [0, 0]
    if n == 1: return [X, Z]
    if n == 2: return double_h(A, B, C, X, Z)
    U, V = X, Z
    T, W = double_h(A, B, C, X, Z)
    bit_length = int(n).bit_length()
    for j in range(bit_length - 2, 0, -1):
        if (n >> j) & 1:
            U, V = add_h(A, B, C, T, W, U, V, X, Z)
            T, W = double_h(A, B, C, T, W)
        else:
            T, W = add_h(A, B, C, U, V, T, W, X, Z)
            U, V = double_h(A, B, C, U, V)
    if n & 1: return add_h(A, B, C, U, V, T, W, X, Z)
    return double_h(A, B, C, U, V)

class ECM(Factorizer):
    def factor(self, B1 = 10000, B2 = 1000000, D = 100):
        assert B1 % 2 == 0 and B2 % 2 == 0
        N = self.N
        while N % 2 == 0:
            N //= 2
        while N % 3 == 0:
            N //= 6
        power_ladder = { p: int(log(B1, p)) for p in get_primes(B1)}
        # Implementation of Algorithm 7.4.4
        def _process_curve(seed): # Here, seed is sigma
            u = ModN(seed*seed - 5, N)
            v = ModN(4*seed, N)
            C = ((v-u)**3 * (3*u+v))/(4*u**3*v) - 2 # Underlying curve is y^2 = x^3 + Cx^2 + x
            # i.e. (A,B,C) = (1,0,C)
            Q_X, Q_Z = u**3, v**3

            # Stage 1
            for p, e in power_ladder.items():
                Q_X, Q_Z = point_multiplication(1,0,C,Q_X,Q_Z, p**e)
            g = gcd(int(Q_Z), N)
            if 1 < g < N: return g
            # Stage 2
            S1_X, S1_Z = double_h(1,0,C,Q_X,Q_Z)
            S2_X, S2_Z = double_h(1,0,C,S1_X,S1_Z)
            S = [None, (S1_X, S1_Z), (S2_X, S2_Z)]
            for d in range(1, D+1):
                if d > 2:
                    Snext_X, Snext_Z = addh(1,0,C,)
