from collections import namedtuple
from typing import Optional

class EllipticCurve(namedtuple("EllpticCurve", ["p", "a", "b"])):
    def __init__(self, p: int, a: int, b: int):
        assert (4*a**3 + 27*b**2) % p != 0
        super().__init__(p,a,b)

class EllipticCurvePoint:
    def __init__(self, curve: EllipticCurve, x: Optional[int], y: Optional[int]):
        p,a,b = curve.p, curve.a, curve.b
        # Point at infinity represented by (x,y) = (None, None)
        assert (x == None and y == None) or (y**2 - (x**3 + a*x + b)) % p == 0
        self.curve = curve
        self.x = x
        self.y = y
    def is_infinity(self):
        return self.x == None and self.y == None
    def __add__(self, other):
        if not(isinstance(other, EllipticCurvePoint)):
            raise ValueError("type error!")
        if self.is_infinity():
            return other
        if other.is_infinity():
            return self
        x1, y1 = self.x, self.y
        x2, y2 = other.x, other.y
        if (x1 == x2 and y2 == -y1)
