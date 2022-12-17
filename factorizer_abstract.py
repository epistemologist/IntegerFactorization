from abc import ABC, abstractmethod
from collections import defaultdict, Counter
from typing import List
from math import prod

from utility import is_prime, Int

class FactorList:
    def __init__(self, l):
        if type(l) == dict:
            self.factor_list = l
        elif type(l) == list:
            self.factor_list = dict(Counter(l))
        else:
            raise ValueError("Invalid data type provided to FactorList!")
    def prod(self):
        return prod({pow(p,e) for p,e in self.factor_list.items()})
    def check_prime(self):
        return all([is_prime(p) for p in self.factor_list])
    def __eq__(self, other):
        return isinstance(other, FactorList) and self.prod() == other.prod()
    def __repr__(self):
        return " . ".join([
            f"{p}" + ("" if e == 1 else f"^{e}")
            for p,e in sorted(self.factor_list.items(),
                              key = lambda x: x[0])
        ])
    def to_dict(self):
        return defaultdict(int, self.factor_list)
class Factorizer(ABC):
    def __init__(self, N):
        self.N = Int(N)

    @abstractmethod
    def factor(self, **args) -> FactorList:
        ...
