import functools
import numpy
import operator
import re
import sympy

from collections import defaultdict
from multiset import FrozenMultiset
from fractions import Fraction as Q
from copy import deepcopy

from .field import Field


class Monomial(FrozenMultiset):
    """A FrozenMultiset representation of a Monomial. Positive integer multiplicities represent powers."""

    def __init__(self, letters_and_powers={}):
        if isinstance(letters_and_powers, (dict, tuple, FrozenMultiset, Monomial)):
            super(Monomial, self).__init__(letters_and_powers)
        elif isinstance(letters_and_powers, str):
            super(Monomial, self).__init__(self.__rstr__(letters_and_powers))
        else:
            print("entry:", repr(letters_and_powers))
            print("type:", type(letters_and_powers))
            raise NotImplementedError

    @staticmethod
    def __rstr__(monomial):
        if monomial == '':
            return dict()
        monomial = monomial.replace("**", "^")
        factors = monomial.split('*')
        factor_groups = defaultdict(list)
        for factor in factors:
            if '^' in factor:
                function, power = factor.split('^')
                power = int(power)
            else:
                function = factor
                power = 1
            factor_groups[function].append(power)
        powers = [(function, sum(powers)) for function, powers in factor_groups.items()]
        return dict(powers)

    def __repr__(self):
        return f"Monomial(\"{str(self)}\")"

    def __str__(self):
        return "*".join([f"{key}^{val}" for key, val in self.items()])

    def tolist(self):
        return [entry for key, val in self.items() for entry in [key, ] * val]

    def subs(self, values_dict):
        return functools.reduce(operator.mul, [values_dict[key] ** val for key, val in self.items()], 1)

    def __mul__(self, other):
        assert isinstance(other, Monomial)
        return super(Monomial, self).__add__(other)

    def __truediv__(self, other):
        raise Exception("Monomial division not implement. Do you mean this to be a Rational Function?")

    def __add__(self, other):
        raise Exception("Monomial addition not implement. Do you mean this to be a Polynomial?")

    def __sub__(self, other):
        raise Exception("Monomial subtraction not implement. Do you mean this to be a Polynomial?")


class Polynomial(object):
    """Generalization of the concept of Multiset where multiplicities are in an arbitrary Field and the elements are Monomials."""

    # This allows the intantiation step to return an element of the udnerlying field.
    # def __new__(cls, coeffs_and_monomials, field):
    #     self = super(Polynomial, cls).__new__(cls)
    #     self.__init__(coeffs_and_monomials, field)
    #     if len(self) == 1 and self.coeffs_and_monomials[0][1] == Monomial(""):
    #         return self.coeffs_and_monomials[0][0]
    #     return self

    def __init__(self, coeffs_and_monomials, field):
        if isinstance(coeffs_and_monomials, str) or isinstance(coeffs_and_monomials, sympy.Basic):
            coeffs_and_monomials = self.__rstr__(str(coeffs_and_monomials), field)
        elif isinstance(coeffs_and_monomials, Polynomial):
            coeffs_and_monomials = coeffs_and_monomials.coeffs_and_monomials
        elif isinstance(coeffs_and_monomials, list):
            assert all([isinstance(entry, tuple) for entry in coeffs_and_monomials])
            assert all([len(entry) == 2 for entry in coeffs_and_monomials])
            assert all([entry[0] in field for entry in coeffs_and_monomials if entry[0] != 0])
            assert all([isinstance(entry[1], Monomial) for entry in coeffs_and_monomials])
        elif coeffs_and_monomials in field:
            coeffs_and_monomials = [(coeffs_and_monomials, Monomial(''))]
        else:
            raise NotImplementedError(f"Received {coeffs_and_monomials} \n of type {type(coeffs_and_monomials)}")
        self.coeffs_and_monomials = coeffs_and_monomials
        self._field = field

    def __str__(self, for_repr=False):
        return " + ".join((f"({str(coeff) if not for_repr else repr(coeff)})"
                           if hasattr(self.field, "name") and self.field.name in ["padic", "finite field", ] else
                           f"{str(coeff) if not for_repr else repr(coeff)}") +
                          (f"*{monomial}" if str(monomial) != "" else "")
                          for coeff, monomial in self.coeffs_and_monomials).replace("+-", "-").replace("+ -", "- ")

    def __repr__(self):
        return f"Polynomial(\"{self.__str__(for_repr=True)}\", {self.field})"

    def __getitem__(self, index):
        return self.coeffs_and_monomials[index]

    @property
    def coeffs_and_monomials(self):
        return self._coeffs_and_monomials

    @coeffs_and_monomials.setter
    def coeffs_and_monomials(self, temp_coeffs_and_monomials):
        self._coeffs_and_monomials = [(coeff, monomial) for coeff, monomial in temp_coeffs_and_monomials if coeff != 0]
        self.reduce()
        if self._coeffs_and_monomials == []:  # ensure there is at least an entry
            self._coeffs_and_monomials = [(0, Monomial(''))]

    @property
    def field(self):
        return self._field

    @field.setter
    def field(self, temp_field):
        if not isinstance(temp_field, Field) and temp_field not in [int, Q]:
            raise Exception(f"temp_field: {temp_field} is not a field")
        if not hasattr(self, "field") or temp_field != self.field:
            self._field = temp_field
            self.coeffs_and_monomials = [(temp_field(coeff), monomial) for coeff, monomial in self.coeffs_and_monomials]

    @property
    def coeffs(self):
        return [coeff for coeff, monomial in self.coeffs_and_monomials]

    @coeffs.setter
    def coeffs(self, temp_coeffs):
        self.coeffs_and_monomials = [(temp_coeff, monomial) for temp_coeff, (coeff, monomial) in zip(temp_coeffs, self.coeffs_and_monomials)]

    @property
    def monomials(self):
        return [monomial for _, monomial in self.coeffs_and_monomials]

    def rationalise(self):
        from pyadic.finite_field import vec_chained_FF_rationalize
        rat_coeffs = vec_chained_FF_rationalize([numpy.array(self.coeffs).astype(int), ], [self.field.characteristic, ]).tolist()
        rat_poly = deepcopy(self)
        rat_poly.coeffs = rat_coeffs
        rat_poly.field = Q
        return rat_poly

    def reduce(self):
        """Merges equal monomials"""
        unequal_monoms = set([monom for _, monom in self.coeffs_and_monomials])
        if len(unequal_monoms) == len(self.coeffs_and_monomials):
            return  # already reduced
        new_coeffs_and_monomials = []
        for unequal_monom in unequal_monoms:
            merged_coeff = sum([coeff for coeff, monom in self.coeffs_and_monomials if monom == unequal_monom])
            new_coeffs_and_monomials += [(merged_coeff, unequal_monom)]
        self.coeffs_and_monomials = new_coeffs_and_monomials

    @staticmethod
    def __rstr__(polynomial, field):
        polynomial = polynomial.replace(" ", "").replace("+-", "-")
        polynomial = re.sub(r"(\+|\-)I\*{0,1}([\d\.]+)", r"\1\2j", polynomial)  # format complex nbrs
        parentheses = [(("(", ), (")", )), (("{", ), ("}", )), (("[", "⟨", "<", ), ("]", "⟩", ">"))]
        lopen_parentheses = [parenthesis[0] for parenthesis in parentheses]
        lclos_parentheses = [parenthesis[1] for parenthesis in parentheses]
        parentheses_balance = [0 for _ in parentheses]
        next_match = ""
        coeffs_and_monomials_strings = []
        for char in polynomial:
            if (char == "+" or char == "-") and all([parenthesis_balance == 0 for parenthesis_balance in parentheses_balance]):
                if next_match != "":
                    coeffs_and_monomials_strings += [next_match.replace(" ", "").replace("**", "^")]
                next_match = char
            else:
                if any([char in open_parentheses for open_parentheses in lopen_parentheses]):
                    parentheses_balance[[char in open_parentheses for open_parentheses in lopen_parentheses].index(True)] += 1
                elif any([char in clos_parentheses for clos_parentheses in lclos_parentheses]):
                    parentheses_balance[[char in clos_parentheses for clos_parentheses in lclos_parentheses].index(True)] -= 1
                next_match += char
        else:
            assert all([parenthesis_balance == 0 for parenthesis_balance in parentheses_balance])
            coeffs_and_monomials_strings += [next_match.replace(" ", "").replace("**", "^")]
        # print(coeffs_and_monomials_strings)
        coeffs_and_monomials = []
        for coeff_and_monomial_string in coeffs_and_monomials_strings:
            coeff = re.findall(r"^[+-]?[\%\+\-\(\)\.j\d/e]*", coeff_and_monomial_string)[0]
            monomial = coeff_and_monomial_string.replace(coeff, "", 1)
            coeff = "+1" if coeff == "+" or coeff == "" else "-1" if coeff == "-" else coeff
            coeff_denom = re.findall(r"/(\d+)", monomial)
            if coeff_denom != []:
                coeff_denom = int(coeff_denom[0])
                monomial = monomial.replace(f"/{coeff_denom}", "")
            else:
                coeff_denom = 1
            if len(monomial) > 0 and monomial[0] == "*":
                monomial = monomial[1:]
            if coeff[:2] == "+(" and coeff[-1:] == ")":
                coeff = coeff[2:-1]
            # print("string:", coeff_and_monomial_string, "\ncoeff", coeff, "\nmonomial:", monomial)
            coeffs_and_monomials += [(field(coeff) / coeff_denom, Monomial(monomial))]
        return coeffs_and_monomials

    def subs(self, base_point, field=None):
        if field is None:
            field = self.field
        base_point = {str(key): val for key, val in base_point.items()}
        new_coeffs_and_monomials = []
        for coeff, monomial in self.coeffs_and_monomials:
            poly = Polynomial(coeff * monomial.subs(base_point), field)
            new_coeffs_and_monomials += [*(poly.coeffs_and_monomials)]
        new_poly = Polynomial(new_coeffs_and_monomials, field)
        new_poly.reduce()
        return new_poly

    def __call__(self, *args, **kwargs):
        return self.subs(*args, **kwargs)

    def __len__(self):
        return len(self.coeffs_and_monomials)

    def __eq__(self, other):
        if other == 0:
            if len(self) == 1 and self.coeffs[0] == 0 and self.monomials[0] == Monomial(''):
                return True
            else:
                return False
        elif isinstance(other, Polynomial):
            return set(self.coeffs_and_monomials) == set(other.coeffs_and_monomials)
        elif isinstance(other, str):
            try:
                return set(self.coeffs_and_monomials) == set(Polynomial(other, self.field).coeffs_and_monomials)
            except:  # noqa
                return False
        else:
            raise NotImplementedError(f"Operation: __eq__; self: {self}; self class {self.__class__}; other: {other}; other class {other.__class__}.")

    def __truediv__(self, other):
        return Polynomial([(coeff / other, monom) for coeff, monom in self.coeffs_and_monomials], self.field)

    def __add__(self, other):
        if isinstance(other, Polynomial):
            return Polynomial(self.coeffs_and_monomials + other.coeffs_and_monomials, self.field)
        else:
            raise NotImplementedError(f"Operation: __add__; self: {self}; self class {self.__class__}; other: {other}; other class {other.__class__}.")

    def __sub__(self, other):
        return self + (-1 * other)

    def __mul__(self, other):
        if isinstance(other, Polynomial):
            new_coeffs_and_monomials = []
            for coeff1, monomial1 in self.coeffs_and_monomials:
                for coeff2, monomial2 in other.coeffs_and_monomials:
                    new_coeffs_and_monomials += [(coeff1 * coeff2, monomial1 * monomial2)]
            return Polynomial(new_coeffs_and_monomials, self.field)
        elif isinstance(other, (int, Q)) or other in self.field:
            new_coeffs_and_monomials = [(other * coeff, monomial) for coeff, monomial in self.coeffs_and_monomials]
            return Polynomial(new_coeffs_and_monomials, self.field)
        else:
            raise NotImplementedError(f"Operation: __mul__; self: {self}; self class {self.__class__}; other: {other}; other class {other.__class__}.")

    def __rmul__(self, other):
        return self * other

    def __neg__(self):
        return -1 * self

    def __pow__(self, n):
        assert (isinstance(n, int) or n.is_integer())
        if n < 0:
            raise Exception("Polynomial to negative power is Rational Function.")
        elif n == 0:
            return Polynomial('1', self.field)
        elif n % 2 == 0:
            root_2_res = self ** (n / 2)
            return root_2_res * root_2_res
        else:
            return self * (self ** (n - 1))
