import sympy

from copy import deepcopy


class RingPoint(dict):
    """Represents a numerical or semi-numericla point in the variety associated to the ring.
    Generalizes the idea of a phase space point from particle physics.
    """

    def __init__(self, ring, field, seed=None):
        super().__init__(ring.random_point(field, seed=seed))
        self.ring, self.field = ring, field

    def __call__(self, string):
        return eval(string, {}, self)

    def univariate_slice(self, seed=None):
        t = sympy.symbols('t')
        self.update(self.ring.univariate_slice(self.field, seed=seed)(t))

    def subs(self, myDict):
        for key, val in self.items():
            self[key] = self.field(val.subs(myDict))

    def copy(self):
        return deepcopy(self)

    def __hash__(self):
        return hash(frozenset(super().items()))
