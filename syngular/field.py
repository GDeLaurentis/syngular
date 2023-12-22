import mpmath
import random

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
        return "({}, {}, {})".format(self.name, self.characteristic, self.digits)

    def __repr__(self):
        return str(self)

    def __call__(self, other):
        """Cast to field."""
        if self.name == "mpc":
            return mpmath.mpc(other)
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

    @property
    def sqrt(self):
        if self.name == "finite field":
            return finite_field_sqrt
        elif self.name == "padic":
            return padic_sqrt
        elif self.name == "mpc":
            return mpmath.sqrt
        else:
            raise Exception(f"Field not understood: {self.field.name}")

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


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #


field = Field()
