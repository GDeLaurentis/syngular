import sympy
import re

from copy import deepcopy

from .ideal import Ideal
from .monomial import non_unicode_powers


class RingPoint(dict):
    """Represents a numerical or semi-numerical point on a variety within the space defiend by the ring.
    Generalizes the idea of a phase space point from particle physics.
    """

    def __init__(self, ring, field, seed=None, val=None, ):
        if val is None:
            super().__init__(ring.random_point(field, seed=seed))
        else:
            super().__init__(val)
        self.ring, self.field = ring, field

    @staticmethod
    def _parse(string):
        string = non_unicode_powers(string)
        string = ' '.join(string.split())
        string = re.sub(r'(\() ([a-zA-Z0-9\-\+])', r'\1\2', string)  # not all white spaces are products, e.g. in '( a'
        string = string.replace(" + ", "+").replace(" - ", "-").replace("^", "**")
        string = string.replace(' ', '*').replace('·', '*').replace(")(", ")*(")
        string = re.sub(r'(\d)([a-zA-Z\(])', r'\1*\2', string)
        string = re.sub(r'([a-zA-Z])(\()', r'\1*\2', string)
        return string

    def __call__(self, string):
        return eval(self._parse(string), {}, self)

    def univariate_slice(self, extra_approximate_constraints=(), indepSet=None, seed=None, verbose=False):
        t = sympy.symbols('t')
        self.update(self.ring.univariate_slice(self.field, extra_approximate_constraints=extra_approximate_constraints,
                                               indepSet=indepSet, seed=seed, verbose=verbose)(t))
        return self

    def subs(self, myDict):
        for key, val in self.items():
            self[key] = self.field(val.subs(myDict))
        return self

    def copy(self):
        return deepcopy(self)

    def __hash__(self):
        return hash((frozenset(super().items()), self.ring, self.field))

    def singular_variety(self, directions_or_ideal=None, valuations=tuple(), seed=None, verbose=False):
        if isinstance(directions_or_ideal, Ideal):
            I = directions_or_ideal
            directions = None
        else:
            I = Ideal(self.ring, directions_or_ideal)
            directions = directions_or_ideal
        self.update(I.point_on_variety(field=self.field, directions=directions, valuations=valuations, seed=seed, verbose=verbose))
        return self
