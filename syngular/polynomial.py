import numbers
import numpy
import re
import sympy
import syngular
import functools
import operator

from collections.abc import Sequence

from fractions import Fraction as Q
from copy import deepcopy

from .field import Field
from .monomial import Monomial


class Polynomial(object):
    """Generalization of the concept of Multiset where multiplicities are in an arbitrary Field and the elements are Monomials."""

    # This allows the intantiation step to return an element of the underlying field.
    # def __new__(cls, coeffs_and_monomials, field):
    #     self = super(Polynomial, cls).__new__(cls)
    #     self.__init__(coeffs_and_monomials, field)
    #     if len(self) == 1 and self.coeffs_and_monomials[0][1] == Monomial(""):
    #         return self.coeffs_and_monomials[0][0]
    #     return self

    def __init__(self, coeffs_and_monomials, field):
        if isinstance(coeffs_and_monomials, (str, int)) or isinstance(coeffs_and_monomials, sympy.Basic):
            coeffs_and_monomials = self.__rstr__(str(coeffs_and_monomials), field)
        elif isinstance(coeffs_and_monomials, Polynomial):
            coeffs_and_monomials = coeffs_and_monomials.coeffs_and_monomials
        elif isinstance(coeffs_and_monomials, (list, tuple)):
            if not all([isinstance(entry, tuple) for entry in coeffs_and_monomials]):
                raise ValueError("Expecting list of tuples, coefficients and monomials.")
            if not all([len(entry) == 2 for entry in coeffs_and_monomials]):
                raise ValueError("Expecting list of tuples, coefficients and monomials, with exactly two entries.")
            if not all([entry[0] in field or entry[0] is None for entry in coeffs_and_monomials if entry[0] != 0]):
                try:
                    coeffs_and_monomials = [(field(coeff), monomial) for coeff, monomial in coeffs_and_monomials]
                except Exception as e:
                    raise ValueError(f"Expecting coefficients to be in the field {field},"
                                     f"received {[coeff for coeff, _  in coeffs_and_monomials]}.\n"
                                     f"Attempted casting to field {field}, got:\n{e}.")
            if not all([isinstance(entry[1], Monomial) for entry in coeffs_and_monomials]):
                raise ValueError("Expecting second entries to be monomials.")
        elif coeffs_and_monomials in field:
            coeffs_and_monomials = [(coeffs_and_monomials, Monomial(''))]
        else:
            raise NotImplementedError(f"Received {coeffs_and_monomials} \n of type {type(coeffs_and_monomials)}")
        self._field = field
        self.coeffs_and_monomials = coeffs_and_monomials

    def __getstate__(self):
        return (tuple(self.coeffs_and_monomials), self.field)

    def __setstate__(self, state):
        self.__init__(*state)

    def __str__(self, for_repr=None):
        if for_repr is None:
            for_repr = not syngular.USE_ELLIPSIS_FOR_PRINT

        def atom_str(coeff_and_monomial, for_repr=False):
            coeff, monomial = coeff_and_monomial
            if coeff is None:
                coeff_bit = "?"
            elif hasattr(self.field, "name") and self.field.name in {"padic", "finite field"}:
                coeff_bit = f"({repr(coeff) if for_repr else str(coeff)})"
            else:
                coeff_bit = str(coeff)
            if not syngular.FORCECDOTS or str(monomial) == '':
                monom_bit = f"{monomial}"
            else:
                monom_bit = f"{syngular.CDOTCHAR}{monomial}"
            return coeff_bit + monom_bit

        terms = self.coeffs_and_monomials
        n = len(terms)

        if for_repr or n <= 3:
            output = "+".join(map(functools.partial(atom_str, for_repr=for_repr), terms))
        else:
            ellipsis = f"...⟪{len(self)-2} terms⟫..."
            output = "+".join([atom_str(self.coeffs_and_monomials[0]), ellipsis, atom_str(self.coeffs_and_monomials[-1])])

        return output.replace("+-", "-")

    def __repr__(self):
        return f"Polynomial(\"{self.__str__(for_repr=True)}\", {self.field})"

    def __getitem__(self, item):
        # Single index or Python slice object
        if isinstance(item, (numbers.Integral, slice)):
            return self.coeffs_and_monomials[item]
        # NumPy boolean mask
        elif (
            (isinstance(item, numpy.ndarray) and item.dtype == bool) or
            (isinstance(item, Sequence) and all(isinstance(i, bool) for i in item))
        ):
            mask = numpy.array(item)
            if mask.shape[0] != len(self):
                raise IndexError(f"Boolean mask length {mask.shape[0]} does not match Polynomial length {len(self)}")
            selected = [x for x, keep in zip(self.coeffs_and_monomials, mask) if keep]
        # NumPy integer array or Python list of integers or list of lists
        elif isinstance(item, (numpy.ndarray, Sequence)) and all(isinstance(i, numbers.Integral) for i in item):
            selected = [self.coeffs_and_monomials[i] for i in item]
        else:
            raise NotImplementedError(f"Unsupported indexing type: {type(item)}")
        # For all non-single-Term results: return a new Terms object with metadata copied
        result = Polynomial(selected, self.field)
        return result

    @property
    def coeffs_and_monomials(self):
        return self._coeffs_and_monomials

    @coeffs_and_monomials.setter
    def coeffs_and_monomials(self, temp_coeffs_and_monomials):
        self._coeffs_and_monomials = [(coeff, monomial) for coeff, monomial in temp_coeffs_and_monomials if coeff != 0]
        self.reduce()
        self._coeffs_and_monomials.sort(key=lambda pair: pair[1])
        if self._coeffs_and_monomials == []:  # ensure there is at least an entry
            self._coeffs_and_monomials = [(self.field(0), Monomial(''))]

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
        return [coeff for coeff, _ in self.coeffs_and_monomials]

    @coeffs.setter
    def coeffs(self, temp_coeffs):
        self.coeffs_and_monomials = [(self.field(temp_coeff), monomial) for temp_coeff, (coeff, monomial) in zip(temp_coeffs, self.coeffs_and_monomials)]

    @property
    def monomials(self):
        return [monomial for _, monomial in self.coeffs_and_monomials]

    @property
    def lead_monomial(self):
        return self.monomials[0]

    @property
    def lead_term(self):
        return Polynomial(self.coeffs_and_monomials[:1], self.field)

    @property
    def variables(self):
        return functools.reduce(operator.or_, [monomial.variables for monomial in self.monomials], set())

    @property
    def linvs(self):
        return [monomial.invs for monomial in self.monomials]

    @property
    def lexps(self):
        return [monomial.exps for monomial in self.monomials]

    def rationalise(self):
        from pyadic.finite_field import vec_chained_FF_rationalize
        rat_coeffs = vec_chained_FF_rationalize([numpy.vectorize(int, otypes='O')(numpy.array(self.coeffs)), ],
                                                [self.field.characteristic ** self.field.digits, ]).tolist()
        rat_poly = deepcopy(self)
        rat_poly._field = Q
        rat_poly.coeffs = rat_coeffs
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
        # print("Polynomial.__rstr__", polynomial)
        polynomial = " ".join(polynomial.split())
        polynomial = polynomial.replace("+-", "-").replace("**", "^").replace(" +", "+").replace("+ ", "+").replace(" -", "-").replace("- ", "-")
        # format complex nbrs - begin
        polynomial = re.sub(r"(\+|\-|^)[iI]\*{0,1}([/\d\.e\+\-]+)", r"\1\2j", polynomial)
        polynomial = re.sub(r"(\+|\-|^)([/\d\.e\+\-]+)\*{0,1}[iI]", r"\1\2j", polynomial)
        # format complex nbrs - end
        parentheses = [(("(", ), (")", )), (("{", ), ("}", )), (("[", "⟨", "<", ), ("]", "⟩", ">"))]
        lopen_parentheses = [parenthesis[0] for parenthesis in parentheses]
        lclos_parentheses = [parenthesis[1] for parenthesis in parentheses]
        parentheses_balance = [0 for _ in parentheses]
        next_match = ""
        coeffs_and_monomials_strings = []
        # print(polynomial)
        for char in polynomial:
            if (char == "+" or char == "-") and all([parenthesis_balance == 0 for parenthesis_balance in parentheses_balance]):
                if next_match != "":
                    coeffs_and_monomials_strings += [next_match]
                next_match = char
            else:
                if any([char in open_parentheses for open_parentheses in lopen_parentheses]):
                    parentheses_balance[[char in open_parentheses for open_parentheses in lopen_parentheses].index(True)] += 1
                elif any([char in clos_parentheses for clos_parentheses in lclos_parentheses]):
                    parentheses_balance[[char in clos_parentheses for clos_parentheses in lclos_parentheses].index(True)] -= 1
                next_match += char
        else:
            if not all([parenthesis_balance == 0 for parenthesis_balance in parentheses_balance]):
                raise AssertionError(f"Unbalanced parentheses: {parentheses_balance}.\n"
                                     f"Match so far: {coeffs_and_monomials_strings}, next one would have been {next_match}")
            coeffs_and_monomials_strings += [next_match]
        # print(coeffs_and_monomials_strings)
        coeffs_and_monomials = []
        for coeff_and_monomial_string in coeffs_and_monomials_strings:
            coeff_pattern = re.compile(r"""
                (?:[?\d\(\)\.\s%/+^-]+        # Common characters
                    | (?<=\d) j               # 'j' if preceeded by digit
                    | e (?=[+-]?\d)           # scientific notation
                    | \* (?=\d)               # '*' if followed by digit
                    | O (?=\(\d)              # O(...) padic term
                    | \^ \d+                  # exponent
                )*
            """, re.VERBOSE)
            # coeff = re.findall(r"^[+-]?[\ \%\+\-\(\)\.j\d/e]*", coeff_and_monomial_string)[0]
            coeff = coeff_pattern.findall(coeff_and_monomial_string)[0]
            while len(coeff) > 0 and coeff[-1] == "(":
                coeff = coeff[:-1]
            monomial = coeff_and_monomial_string.replace(coeff, "", 1)
            # print(f"coeff: {coeff}\nmonom: {monomial}")
            coeff = coeff.replace(" ", "")
            coeff = "+1" if coeff == "+" or coeff == "" else "-1" if coeff == "-" else coeff
            coeff_denom = re.findall(r"/(\d+)$", monomial)
            if coeff_denom != []:
                coeff_denom = int(coeff_denom[0])
                monomial = monomial.replace(f"/{coeff_denom}", "")
            else:
                coeff_denom = 1
            if len(monomial) > 0 and monomial[0] == "*":
                monomial = monomial[1:]
            if coeff[:2] == "+(" and coeff[-1:] == ")":
                coeff = coeff[2:-1]
            if coeff[:1] == "(" and coeff[-1:] == ")":
                coeff = coeff[1:-1]
            # print("string:", coeff_and_monomial_string, "\ncoeff:", coeff, "\ncoeff denom:", coeff_denom, "\nmonomial:", monomial)
            coeffs_and_monomials += [(field(coeff) / coeff_denom if coeff not in ('?', '+?', '-?') else None, Monomial(monomial))]
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
        if isinstance(other, Monomial):
            return Polynomial([(coeff, monom / other) for coeff, monom in self.coeffs_and_monomials], self.field)
        if isinstance(other, Polynomial):
            if len(self) == len(other) == 1 and other.monomials[0].issubset(self.monomials[0]):  # trivial div e.g. between lead_terms
                return Polynomial([(self_coeff / other_coeff, self_monom / other_monom) for (self_coeff, self_monom), (other_coeff, other_monom)
                                   in zip(self.coeffs_and_monomials, other.coeffs_and_monomials)], self.field)
            else:
                quotient, remainder = divmod(self, other)
                if remainder == 0:
                    return quotient
                return NotImplemented
        return Polynomial([(coeff / other, monom) for coeff, monom in self.coeffs_and_monomials], self.field)

    def __divmod__(self, other):
        if not isinstance(other, Polynomial):
            return NotImplemented

        dividend = deepcopy(self)
        divisor = deepcopy(other)
        quotients = []

        while dividend != 0:
            if not divisor.lead_monomial.issubset(dividend.lead_monomial):
                break  # divisor does not divide dividend, stop
            quotient = dividend.lead_term / divisor.lead_term
            quotients.append(quotient)
            previous_dividend = deepcopy(dividend)
            dividend -= quotient * divisor
            if dividend == previous_dividend:
                raise Exception(f"Failed to compute __divmod__ between\n{self}\nand\n{other}\nPerhaps coefficients are O(p^ν)?")

        quotient = sum(quotients, Polynomial('0', self.field))
        remainder = dividend

        return quotient, remainder

    def __floordiv__(self, other):
        quotient, _ = divmod(self, other)
        return quotient

    def __mod__(self, other):
        _, remainder = divmod(self, other)
        return remainder

    def __add__(self, other):
        if isinstance(other, Polynomial):
            return Polynomial(self.coeffs_and_monomials + other.coeffs_and_monomials, self.field)
        elif isinstance(other, self.field.random().__class__):
            return Polynomial(self.coeffs_and_monomials + [(other, Monomial('')), ], self.field)
        else:
            try:
                return Polynomial(self.coeffs_and_monomials + [(self.field(other), Monomial('')), ], self.field)
            except ValueError:
                return NotImplemented
            # for debugging
            # raise NotImplementedError(f"Operation: __add__; self: {self}; self class {self.__class__}; other: {other}; other class {other.__class__}.")

    def __radd__(self, other):
        return self + other

    def __sub__(self, other):
        return self + (-1 * other)

    def __mul__(self, other):
        if isinstance(other, Monomial):
            other = Polynomial([(1, other)], self.field)
        if isinstance(other, Polynomial):
            new_coeffs_and_monomials = []
            for coeff1, monomial1 in self.coeffs_and_monomials:
                for coeff2, monomial2 in other.coeffs_and_monomials:
                    new_coeffs_and_monomials += [(None if coeff1 is None or coeff2 is None else coeff1 * coeff2, monomial1 * monomial2)]
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
