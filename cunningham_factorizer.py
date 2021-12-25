from collections import namedtuple

from factorizer_abstract import FactorList

class CunninghamNumber(namedtuple("CunninghamNumber", ["b", "n", "k"])):
    # Class to represent numbers of the form b^n +/- 1
    def __init__(self, b, n, k):
        assert abs(k) == 1
        super().__init__(b,n,k)
    def __int__(self):
        return pow(self.b, self.n) + self.k
    def __repr__(self):
        sign = "+" if self.k > 0 else "-"
        return f"{self.b}^{self.n}" + (f" {sign} {abs(self.k)}" if self.k != 0 else "")


