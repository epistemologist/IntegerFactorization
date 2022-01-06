from collections import namedtuple
from typing import Optional, Tuple
from random import randint
from utility import ModN, GMPY_IMPORT, tonelli_shanks, legendre

class EllipticCurve(namedtuple("EllpticCurve", ["p", "a", "b"])):
    def __init__(self, p: int, a: int, b: int):
        # Assert curve is not singular
        assert (4*a**3 + 27*b**2) % p != 0
        super().__init__()
    def random_point(self) -> Tuple[int, int]:
        while True:
            x = randint(0, self.p-1)
            t = (x * (x*x + self.a) + self.b) % self.p
            if legendre(t, self.p) == -1: continue
            return (x, tonelli_shanks(t, self.p))
    def get_point(self, x: int) -> Optional[Tuple[int, int]]:
        t = (x * (x*x + self.a) + self.b) % self.p
        if legendre(t, self.p) == -1: return None
        return (x, tonelli_shanks(t, self.p))

class EllipticCurveAffinePoint:
    def __init__(self, curve: EllipticCurve, x,y, is_infinite=False):
        self.curve = curve
        self.x = ModN(x, curve.p)
        self.y = ModN(y, curve.p)
        self.is_infinite = is_infinite
        a,b = curve.a, curve.b
        assert self.is_infinite or self.x**3 + a*self.x + b == self.y**2
    def __repr__(self):
        return "O" if self.is_infinite else f"({self.x}, {self.y})"
    def __iter__(self):
        return iter((self.x, self.y, self.is_infinite))
    def __add__(self, other):
        assert type(other) == EllipticCurveAffinePoint
        return affine_add(self, other)
    def __sub__(self, other):
        assert type(other) == EllipticCurveAffinePoint
        return affine_add(self, affine_neg(other))
    def __mul__(self, other):
        assert type(other) == int
        return affine_mul(self, other)
    def __rmul__(self, other):
        assert type(other) == int
        return affine_mul(self, other)

# See "Prime Numbers: A Computational Perspective", Algorithm 7.2.2 for implementation details

def affine_add(p1: EllipticCurveAffinePoint, p2: EllipticCurveAffinePoint) -> EllipticCurveAffinePoint:
    x1, y1, p1_infinite = tuple(p1)
    x2, y2, p2_infinite = tuple(p2)
    assert p1.curve == p2.curve
    a = p1.curve.a
    if p1_infinite: return p2
    if p2_infinite: return p1
    if x1 == x2:
        if y1 + y2 == 0: return EllipticCurveAffinePoint(0,0,is_infinite=True)
        m = (3*x1*x1 + a) / (2*y1)
    else:
        m = (y2 - y1) / (x2 - x1)
    x3 = m*m - x1 - x2
    return EllipticCurveAffinePoint(p1.curve, x3, m*(x1-x3)-y1)

def affine_neg(p: EllipticCurveAffinePoint) -> EllipticCurveAffinePoint:
    x, y, is_infinite = tuple(p)
    return EllipticCurvePoint(p.curve,x,-y,is_infinite=is_infinite)

def affine_mul(x: EllipticCurveAffinePoint, n: int) -> EllipticCurveAffinePoint:
    # Multiplication by doubling
    y = EllipticCurveAffinePoint(x.curve, 0,0,is_infinite=True)
    if n == 0: return y
    while n > 1:
        if n % 2 == 0:
            x = x+x
            n //= 2
        else:
            y = x+y
            x = x+x
            n = (n-1) // 2
    return x+y


