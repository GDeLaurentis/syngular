import mpmath
import random

from fractions import Fraction

from pyadic import PAdic, ModP, GaussianRational
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

    def __getstate__(self):
        return (self.name, self.characteristic, self.digits)

    def __setstate__(self, state):
        self.set(*state)

    def __eq__(self, other):
        return isinstance(other, Field) and self.name == other.name and self.characteristic == other.characteristic and self.digits == other.digits

    def __str__(self):
        return f"Field('{self.name}', {self.characteristic}, {self.digits})"

    __repr__ = __str__

    def __call__(self, *args):
        """Cast to field."""
        if self.name in ["mpc", "C"]:
            return mpmath.mpc(mpmath.mpmathify(args[0]))
        elif self.name in ["padic", "Qp"]:
            return PAdic(args[0], self.characteristic, self.digits)
        elif self.name in ["finite field", "Fp"]:
            return ModP(args[0], self.characteristic)
        elif self.name in ["rational", "Q"]:
            return Fraction(*args)
        elif self.name in ["gaussian rational", "Q[i]", "Qi"]:
            return GaussianRational(*args)
        else:
            raise NotImplementedError

    @property
    def name(self):
        return self._name

    @name.setter
    def name(self, value):
        if value not in ['rational', 'Q', 'mpc', 'C', 'gaussian rational', 'Q[i]', 'Qi', 'finite field', 'Fp', 'padic', 'Qp']:
            raise Exception("Field must be one of 'mpc', 'gaussian rational', 'finite field', 'padic', or aliases thereof.")
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
        if self.name in ['mpc', 'C']:
            return mpmath.mp.dps
        else:
            return self._digits

    @digits.setter
    def digits(self, value):
        if value < 0:
            raise Exception("Digits must be positive.")
        elif self.name in ['mpc', 'C']:
            mpmath.mp.dps = value
        else:
            self._digits = value

    def __contains__(self, other):
        return isinstance(other, self(1).__class__)

    @property
    def is_algebraically_closed(self):
        if self.name in ['mpc', 'C']:
            return True
        else:
            return False

    @property
    def tollerance(self):
        if self.name in ['gaussian rational', 'finite field']:
            return 0
        elif self.name in ['mpc', 'C']:
            return mpmath.mpf('10e-{}'.format(int(min([0.95 * mpmath.mp.dps, mpmath.mp.dps - 4]))))
        elif self.name in ['padic', 'Qp']:
            return PAdic(0, self.characteristic, 0, self.digits)

    @property
    def singular_notation(self):
        if self.name in ['mpc', 'C']:
            return '(complex,{},I)'.format(self.digits - 5)
        elif self.name in ['finite field', 'Fp', 'padic', 'Qp']:
            return str(self.characteristic)
        else:
            return None

    def sqrt(self, val):
        if self.name in ["finite field", 'Fp']:
            if not isinstance(val, ModP):
                val = self(val)
            return finite_field_sqrt(val)
        elif self.name in ["padic", 'Qp']:
            if not isinstance(val, PAdic):
                val = self(val)
            return padic_sqrt(val)
        elif self.name in ["mpc", 'C']:
            return mpmath.sqrt(val)
        else:
            raise Exception(f"Field not understood: {self.field.name}")

    @property
    def one(self):
        return self(1)

    @property
    def zero(self):
        return self(0)

    @property
    def j(self):
        return self.sqrt(-1)

    i = I = j

    def random(self, shape=(1, )):
        if shape == (1, ):
            if self.name in ["padic", 'Qp']:
                p, k = self.characteristic, self.digits
                return PAdic(random.randrange(0, p ** k - 1), p, k)
            elif self.name in ["finite field", 'Fp']:
                p = self.characteristic
                return ModP(random.randrange(0, p), p)
            elif self.name in ["mpc", 'C']:
                return mpmath.mpc(str(Fraction(random.randrange(-100, 101), random.randrange(1, 201))),
                                  str(Fraction(random.randrange(-100, 101), random.randrange(1, 201))))
            else:
                raise NotImplementedError
        else:
            raise NotImplementedError

    def random_square(self):
        from pyadic.field_extension import FieldExtension
        for i in range(10000):
            val = self.random()
            if not isinstance(self.sqrt(val), FieldExtension):
                return val

    def random_element(self, *args, **kwargs):
        import warnings
        warnings.warn(
            "random_element is deprecated and will be removed in a future version. "
            "Use random instead.",
            DeprecationWarning,
            stacklevel=2
        )
        return self.random(*args, **kwargs)

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
            self._ε = self.epsilon()
        return self._ε

    @ε.setter
    def ε(self, temp_ε):
        self._ε = temp_ε


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #


Q = Field("rational", 0, 0)
Qi = Field("gaussian rational", 0, 0)
