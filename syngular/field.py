import mpmath
import random
import re

from fractions import Fraction as Q

from pyadic import PAdic, ModP
from pyadic.padic import padic_sqrt
from pyadic.finite_field import finite_field_sqrt

mpmath.mp.dps = 300


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #


class Field(object):
    """A class representing number fields."""

    def __init__(self, *args):
        """Defaults to ('mpc', 0, 300)"""
        if args == ():
            args = ('mpc', 0, 300)
        self.set(*args)

    def set(self, *args):
        """(name, characteristic, digits)"""
        self.name = args[0]
        self.characteristic = args[1]
        self.digits = args[2]

    def __str__(self):
        return "('{}', {}, {})".format(self.name, self.characteristic, self.digits)

    def __repr__(self):
        return str(self)

    def __call__(self, other):
        """Cast to field."""
        if self.name == "mpc":
            return mpmath.mpc(mpmath.mpmathify(other))
        elif self.name == "padic":
            return PAdic(other, self.characteristic, self.digits)
        elif self.name == "finite field":
            return ModP(other, self.characteristic)
        else:
            raise NotImplementedError

    @property
    def name(self):
        return self._name

    @name.setter
    def name(self, value):
        if value not in ['mpc', 'gaussian rational', 'finite field', 'padic']:
            raise Exception("Field must be one of 'mpc', 'gaussian rational', 'finite field', 'padic'.")
        else:
            self._name = value

    @property
    def characteristic(self):
        return self._characteristic

    @characteristic.setter
    def characteristic(self, value):
        if value < 0:
            raise Exception("Characteristic must be non-negative.")
        else:
            self._characteristic = value

    @property
    def digits(self):
        if self.name == 'mpc':
            return mpmath.mp.dps
        else:
            return self._digits

    @digits.setter
    def digits(self, value):
        if value < 0:
            raise Exception("Digits must be positive.")
        elif self.name == 'mpc':
            mpmath.mp.dps = value
        else:
            self._digits = value

    @property
    def is_algebraically_closed(self):
        if self.name == 'mpc':
            return True
        else:
            return False

    @property
    def tollerance(self):
        if self.name in ['gaussian rational', 'finite field']:
            return 0
        elif self.name == 'mpc':
            return mpmath.mpf('10e-{}'.format(int(min([0.95 * mpmath.mp.dps, mpmath.mp.dps - 4]))))
        elif self.name == 'padic':
            return PAdic(0, self.characteristic, 0, self.digits)

    @property
    def singular_notation(self):
        if self.name == 'mpc':
            return '(complex,{},I)'.format(self.digits - 5)
        elif self.name in ['finite field', 'padic']:
            return str(self.characteristic)
        else:
            return None

    def sqrt(self, val):
        if self.name == "finite field":
            if not isinstance(val, ModP):
                val = self(val)
            return finite_field_sqrt(val)
        elif self.name == "padic":
            if not isinstance(val, PAdic):
                val = self(val)
            return padic_sqrt(val)
        elif self.name == "mpc":
            return mpmath.sqrt(val)
        else:
            raise Exception(f"Field not understood: {self.field.name}")

    @property
    def j(self):
        return self.sqrt(-1)

    i = I = j

    def random_element(self, shape=(1, )):
        if shape == (1, ):
            if self.name == "padic":
                p, k = self.characteristic, self.digits
                return PAdic(random.randrange(0, p ** k - 1), p, k)
            elif self.name == "finite field":
                p = self.characteristic
                return ModP(random.randrange(0, p), p)
            elif self.name == "mpc":
                return mpmath.mpc(str(Q(random.randrange(-100, 101), random.randrange(1, 201))),
                                  str(Q(random.randrange(-100, 101), random.randrange(1, 201))))
            else:
                raise NotImplementedError
        else:
            raise NotImplementedError

    def epsilon(self, shape=(1, ), ):
        if shape == (1, ):
            if not hasattr(self, "_ε"):
                if self.name == "padic":
                    self._ε = self.characteristic
                elif self.name == "finite field":
                    raise ValueError("Finite field infinitesimal does not exist.")
                elif self.name == "mpc":
                    self._ε = mpmath.mpf('1e-30')
                else:
                    raise NotImplementedError
            return self._ε
        else:
            raise NotImplementedError

    @property
    def ε(self):
        if not hasattr(self, "_ε"):
            self._ε = self.epsilon
        return self._ε

    @ε.setter
    def ε(self, temp_ε):
        self._ε = temp_ε


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #


field = Field()
