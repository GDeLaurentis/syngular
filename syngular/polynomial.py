import functools
import numpy
import operator
import re

from collections import defaultdict
from multiset import FrozenMultiset
from fractions import Fraction as Q
from copy import deepcopy

from .field import Field


class Monomial(FrozenMultiset):

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

    def eval(self, values_dict):
        return functools.reduce(operator.mul, [values_dict[key] ** val for key, val in self.items()], 1)


class Polynomial(object):
    
    def __init__(self, coeffs_and_monomials, field):
        if isinstance(coeffs_and_monomials, str):
            coeffs_and_monomials = self.__rstr__(coeffs_and_monomials, field)
        self.coeffs_and_monomials = coeffs_and_monomials
        self._field = field

    def __str__(self, for_repr=False):
        return " + ".join((f"({str(coeff) if not for_repr else repr(coeff)})" if hasattr(self.field, "name") and self.field.name == "padic" else 
                           f"{str(coeff) if not for_repr else repr(coeff)}") + 
                          (f"*{monomial}" if str(monomial) != "" else "") 
                          for coeff, monomial in self.coeffs_and_monomials).replace("+-", "-")
    
    def __repr__(self):
        return f"Polynomial(\"{self.__str__(for_repr=True)}\", {self.field})"

    @property
    def coeffs_and_monomials(self):
        return self._coeffs_and_monomials

    @coeffs_and_monomials.setter
    def coeffs_and_monomials(self, temp_coeffs_and_monomials):
        self._coeffs_and_monomials = [(coeff, monomial) for coeff, monomial in temp_coeffs_and_monomials if coeff != 0]

    @property
    def field(self):
        return self._field

    @field.setter
    def field(self, temp_field):
        if not isinstance(temp_field, Field) and not temp_field in [int, Q]:
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

    def rationalise(self):
        from pyadic.finite_field import vec_chained_FF_rationalize
        rat_coeffs = vec_chained_FF_rationalize([numpy.array(self.coeffs).astype(int), ], [self.field.characteristic, ]).tolist()
        rat_poly = deepcopy(self)
        rat_poly.coeffs = rat_coeffs
        rat_poly.field = Q
        return rat_poly

    def reduce(self):
        """Merges equal monomials"""
        unequal_monoms = set([monom for coeff, monom in self.coeffs_and_monomials])
        new_coeffs_and_monomials = []
        for unequal_monom in unequal_monoms:
            merged_coeff = sum([coeff for coeff, monom in self.coeffs_and_monomials if monom == unequal_monom])
            new_coeffs_and_monomials += [(merged_coeff, unequal_monom)]
        self.coeffs_and_monomials = new_coeffs_and_monomials

    @staticmethod
    def __rstr__(polynomial, field):
        parentheses = [(("(", ), (")", )), (("{", ), ("}", )), (("[", "⟨", "<", ), ("]", "⟩", ">"))]
        lopen_parentheses = [parenthesis[0] for parenthesis in parentheses]
        lclos_parentheses = [parenthesis[1] for parenthesis in parentheses]
        parentheses_balance = [0 for parenthesis in parentheses]
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
            coeff = re.findall(r"^[+-]?[\+\-\(\)\.j\d/]*", coeff_and_monomial_string)[0]
            monomial = coeff_and_monomial_string.replace(coeff, "", 1)
            coeff = "+1" if coeff == "+" or coeff == "" else "-1" if coeff == "-" else coeff
            if len(monomial) > 0 and monomial[0] == "*":
                monomial = monomial[1:]
            if coeff[:2] == "+(" and coeff[-1:] == ")":
                coeff = coeff[2:-1]
            # print("string:", coeff_and_monomial_string, "\ncoeff", coeff, "\nmonomial:", monomial)
            coeffs_and_monomials += [(field(coeff), Monomial(monomial))]
        return coeffs_and_monomials

    def subs(self, base_point, field):
        indep_variables = {str(key): val for key, val in base_point.items() if key != val}
        dep_variables = {str(key): val for key, val in base_point.items() if key == val}
        new_coeffs_and_monomials = []
        for coeff, monomial in self.coeffs_and_monomials:
            dep_part = {key: val for key, val in monomial.items() if key in map(str, dep_variables)}
            indep_part = {key: val for key, val in monomial.items() if key in map(str, indep_variables)}
            new_monomial = Monomial(dep_part)
            new_coeff = coeff * Monomial(indep_part).eval(indep_variables)
            new_coeffs_and_monomials += [(new_coeff, new_monomial)]
        new_poly = Polynomial(new_coeffs_and_monomials, field)
        new_poly.reduce()
        return new_poly

    def __truediv__(self, other):
        return Polynomial([(coeff / other, monom) for coeff, monom in self.coeffs_and_monomials], self.field)

    def __eq__(self, other):
        if other == 0:
            if self.coeffs_and_monomials == []:
                return True
            else:
                return False
        elif isinstance(other, Polynomial):
            return set(self.coeffs_and_monomials) == set(other.coeffs_and_monomials)
        else:
            raise NotImplementedError

    def __len__(self):
        return len(self.coeffs_and_monomials)
