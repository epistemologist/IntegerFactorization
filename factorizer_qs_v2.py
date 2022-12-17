from factorizer_abstract import Factorizer
from dataclasses import dataclass
from utility import Int, ModN, tonelli_shanks
from typing import Tuple

@dataclass
class Polynomial:
    a: Int
    b: Int
    c: Int

    def eval(self, x: Int) -> Int:
        return a*x*x + b*x + c

    def get_mod_roots(self, p: Int) -> Tuple[Int, Int]:
        a, b, c =

class QuadraticSieve(Factorizer):
    def factor(self):
