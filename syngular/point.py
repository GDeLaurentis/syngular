import sympy

from copy import deepcopy

from .ideal import Ideal


class RingPoint(dict):
    """Represents a numerical or semi-numericla point on a variety within the space defiend by the ring.
    Generalizes the idea of a phase space point from particle physics.
    """

    def __init__(self, ring, field, seed=None):
        super().__init__(ring.random_point(field, seed=seed))
        self.ring, self.field = ring, field

    def __call__(self, string):
        return eval(string, {}, self)

    def univariate_slice(self, indepSet=None, seed=None, verbose=False):
        t = sympy.symbols('t')
        self.update(self.ring.univariate_slice(self.field, indepSet=indepSet, seed=seed, verbose=verbose)(t))

    def subs(self, myDict):
        for key, val in self.items():
            self[key] = self.field(val.subs(myDict))

    def copy(self):
        return deepcopy(self)

    def __hash__(self):
        return hash(frozenset(super().items()))

    def singular_variety(self, directions=None, valuations=tuple()):  # to be improved
        I = Ideal(self.ring, directions)
        self.update(I.point_on_variety(field=self.field, directions=directions, valuations=valuations))
