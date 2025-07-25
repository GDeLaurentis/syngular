import functools
import re
import sympy
import numpy
import inspect
import random
import warnings

from packaging.version import Version

from .tools import execute_singular_command, Singular_version
from .ring import Ring
from .qring import QuotientRing
from .ideal_algorithms import Ideal_Algorithms
from .variety import Variety_of_Ideal
from .polynomial import Polynomial
from .field import Field


class Ideal(Ideal_Algorithms, Variety_of_Ideal, object):

    def __init__(self, ring, generators):
        if isinstance(generators, tuple):
            generators = list(generators)
        if not isinstance(ring, Ring) or not isinstance(generators, list):
            raise Exception("Invalid Ideal initialisation.")
        self.ring = ring
        self.generators = generators
        self.test_valid_ideal()
        self._dim = None
        self._indepSets = None

    def test_valid_ideal(self):
        if any(isinstance(entry, (list, tuple)) for entry in self.generators):
            raise ValueError("Ideal generators should not contain lists or tuples. Have you forgotten to flatten?")
        singular_commands = [f"ring r = {self.ring};",
                             f"ideal i = {self};",
                             "print(i);"
                             "$"]
        execute_singular_command(singular_commands)

    @property
    def generators(self):
        return self._generators

    @generators.setter
    def generators(self, generators):
        if not isinstance(generators, list):
            raise Exception("Ideal generators must be list")
        if len(generators) == 0:
            generators = ['0']
        self._generators = generators

    def __hash__(self):
        return hash(", ".join(self.reduced_groebner_basis)) + hash(self.ring)

    def __eq__(self, other):
        return self.reduced_groebner_basis == other.reduced_groebner_basis

    def squash(self):
        self.generators = self.minbase

    def squoosh(self):
        warnings.warn(
            "The 'squoosh' method is deprecated, please use 'squash' instead.",
            DeprecationWarning,
            stacklevel=2
        )
        return self.squash()

    @property
    def dim(self):
        if self._dim is None:
            singular_commands = [f"ring r = {self.ring};",
                                 f"ideal gb = {','.join(self.groebner_basis)};",
                                 # f"ideal i = {self};",   # Check which is better
                                 # "ideal gb = std(i);",   # {','.join(self.groebner_basis)};",
                                 "print(dim(gb));",
                                 "$"]
            output = execute_singular_command(singular_commands)
            self._dim = int(output)
        return self._dim

    @dim.setter
    def dim(self, val):
        self._dim = val

    @property
    def dims(self):
        return {entry.count(1) for entry in self.indepSets}

    @property
    def codim(self):
        return Ideal(self.ring, ['0']).dim - self.dim

    @property
    def codims(self):
        zeroIdeal_dim = Ideal(self.ring, ['0']).dim
        return {zeroIdeal_dim - dim for dim in self.dims}

    def guess_indep_set(self):
        equations = self.generators
        if isinstance(self.ring, QuotientRing):
            equations += self.ring.ideal.generators
        n = len(self.ring.variables)
        m = self.dim if self._dim is not None else (n - len(equations))
        lst = [0] * (n - m) + [1] * m
        equations_variables = [Polynomial(equation, field=Field("rational", 0, 0)).variables for equation in equations]
        for _ in range(1000):  # discard obviously wrong indep sets: all equations must have at least 1 dependent variable
            random.shuffle(lst)
            guess_dict = dict(zip(map(str, self.ring.variables), lst))
            if all([0 in [guess_dict[variable] for variable in variables] for variables in equations_variables]):
                break
        return tuple(lst)

    @functools.cached_property
    def indepSet(self):
        singular_commands = [f"ring r = {self.ring};",
                             f"ideal gb = {','.join(self.groebner_basis)};",
                             "print(indepSet(gb));",
                             "$"]
        output = execute_singular_command(singular_commands)
        return tuple(map(int, output.split(",\n")))

    @property
    def indepSets(self):
        if self._indepSets is None:
            singular_commands = [f"ring r = {self.ring};",
                                 # f"ideal gb = {','.join(self.groebner_basis)};",   # this breaks singular variety construction, especially with mpcs
                                 f"ideal i = {self};",
                                 "ideal gb = groebner(i);",
                                 "print(indepSet(gb, 1));",
                                 "$"]
            output = execute_singular_command(singular_commands)
            if output == 'empty list':
                return [self.indepSet]
            indepSets = [tuple(map(int, line.replace(" ", "").split(","))) for line in output.split("\n") if ":" not in line]
            self._indepSets = indepSets
        return self._indepSets

    @indepSets.setter
    def indepSets(self, val):
        self._indepSets = val

    def get_groebner_basis(self, reduced=False, algorithm=['groebner', 'slimgb'][1]):
        singular_commands = [f"ring r = {self.ring};",
                             f"ideal i = {self};",
                             f"ideal gb = {algorithm}(i);",
                             "short=0;",
                             "print(gb);",
                             "$"]
        if reduced:
            singular_commands = ["option(redSB);"] + singular_commands
        singular_command = "\n".join(singular_commands)
        output = execute_singular_command(singular_command)
        output = [line.replace(",", "") for line in output.split("\n")]
        return output

    @functools.cached_property
    def groebner_basis(self):
        return self.get_groebner_basis(reduced=True, algorithm='groebner')

    @functools.cached_property
    def reduced_groebner_basis(self):
        return self.get_groebner_basis(reduced=True, algorithm='groebner')

    @functools.cached_property
    def leadGBmonomials(self):
        """Gives the leading monomials of the Groebner basis polynomials."""
        singular_commands = [f"ring r = {self.ring};",
                             f"ideal gb = {','.join(self.groebner_basis)};",
                             "def lts = lead(gb);",
                             "short=0;",
                             "print(lts);"]
        output = execute_singular_command(singular_commands)
        # print(output)
        output = [line.replace(",", "") for line in output.split("\n")]
        return output

    @functools.cached_property
    def minbase(self):
        singular_commands = [f"ring r = {self.ring};",
                             f"ideal i = {self};",
                             "short=0;",
                             "print(minbase(i));",
                             "$"]
        output = execute_singular_command(singular_commands)
        output = [line.replace(",", "") for line in output.split("\n")]
        return output

    @functools.cached_property
    def radical(self):
        """Returns the radical of the ideal."""
        singular_commands = ["LIB \"primdec.lib\";",
                             f"ring r = {self.ring};",
                             f"ideal i = {self};",
                             "def pr = radical(i);",   # options: GTZ / SY
                             "short=0;",
                             "print(pr);",
                             "$"]
        output = execute_singular_command(singular_commands)
        output = [line.replace(",", "") for line in output.split("\n")]
        cls, ring = self.__class__, self.ring
        return cls(ring, output)

    @functools.cached_property
    def primary_decomposition(self):
        singular_commands = ["LIB \"primdec.lib\";",
                             f"ring r = {self.ring};",
                             f"ideal i = {self};",
                             # f"ideal gb = {','.join(self.groebner_basis)};",
                             "def pr = primdecGTZ(i);",   # options: GTZ / SY
                             "short=0;",
                             "print(pr);",
                             "$"]
        output = execute_singular_command(singular_commands)
        output = [line.replace(",", "") for line in output.split("\n")]

        def clean_up(string):
            string = re.sub(r"_\[\d+\]=", "", string)
            string = re.sub(r"\[\d+\]:", "|", string)
            return string.replace(" ", "")

        primary_decomposed = ",".join(list(map(clean_up, output))).replace("|,", "|").replace(",|", "|").replace("|,", "|").split("||")
        primary_decomposed = list(filter(lambda x: x != '', primary_decomposed))
        primary_decomposed = [entry.split("|") for entry in primary_decomposed]
        primary_decomposed = [(entry[0].split(","), entry[1].split(",")) for entry in primary_decomposed]

        cls, ring = self.__class__, self.ring
        for i, (primary, prime) in enumerate(primary_decomposed):
            primary_decomposed[i] = (cls(ring, primary), cls(ring, prime))

        return primary_decomposed

    def eliminate(self, var_range):
        singular_commands = [f"ring r1 = {self.ring};",
                             f"ideal i = {self};",
                             f"ideal j = eliminate(i, {var_range.start} .. {var_range.stop});",
                             "print(minbase(j));",
                             "$"]
        output = execute_singular_command(singular_commands)
        output = [line.replace(",", "") for line in output.split("\n")]
        cls, ring = self.__class__, self.ring
        ring = Ring(ring.field, ring.variables[:var_range.start] + ring.variables[var_range.stop:], ring.ordering if isinstance(ring.ordering, str) else 'dp')
        return cls(ring, output)

    def __contains__(self, other):
        """Implements ideal membership."""
        if isinstance(other, Ideal):
            return self.__ideal_contains__(other)
        elif isinstance(other, (str, sympy.Expr, Polynomial)):
            return self.__poly_contains__(other)
        else:
            raise NotImplementedError

    def __ideal_contains__(self, other):
        qIdeal = self / other
        if list(map(str, qIdeal.generators)) == ['1']:
            return True
        else:
            return False

    def __poly_contains__(self, other):
        if not isinstance(other, str):
            other = str(other)
        assert isinstance(other, str)
        singular_commands = [f"ring r = {self.ring};",
                             f"ideal gb = {','.join(self.groebner_basis)};",
                             "poly f = " + other + ";",
                             "print(reduce(f, gb));"
                             "$"]
        output = execute_singular_command(singular_commands)
        output = [line.replace(",", "") for line in output.split("\n")]
        if output == ['0']:
            return True
        else:
            return False

    def __add__(self, other):
        """Addition of Ideals = Intersection of Varieties."""
        singular_commands = [f"ring r = {self.ring};",
                             f"ideal i = {self};",
                             f"ideal j = {other};",
                             "ideal k = i + j;",
                             "print(k);",
                             "$"]
        output = execute_singular_command(singular_commands)
        output = [line.replace(",", "") for line in output.split("\n")]
        cls, ring = self.__class__, self.ring
        return cls(ring, output)

    def __mul__(self, other):
        """Product of Ideals"""
        singular_commands = [f"ring r = {self.ring};",
                             f"ideal i = {self};",
                             f"ideal j = {other};",
                             "ideal k = i * j;",
                             "print(k);",
                             "$"]
        output = execute_singular_command(singular_commands)
        output = [line.replace(",", "") for line in output.split("\n")]
        cls, ring = self.__class__, self.ring
        return cls(ring, output)

    def __pow__(self, n):
        assert type(n) is int and n >= 0
        cls, ring = self.__class__, self.ring
        if n == 0:
            return cls(ring, ['1'])
        elif n % 2 == 0:
            root_2_res = self ** (n // 2)
            return root_2_res * root_2_res
        else:
            return self * (self ** (n - 1))

    def __and__(self, other):
        """Intersection of Ideals = Union of Varieties - This uses Python's set intersection operator '&'."""
        singular_commands = [f"ring r = {self.ring};",
                             f"ideal i = {self};",
                             f"ideal j = {other};",
                             "ideal k = intersect(i, j);",
                             "print(k);",
                             "$"]
        output = execute_singular_command(singular_commands)
        output = [line.replace(",", "") for line in output.split("\n")]
        cls, ring = self.__class__, self.ring
        return cls(ring, output)

    @staticmethod
    def intersection(*args):
        """Intersection of Ideals - wrapper around & operator for chained intersection."""
        if len(args) == 1:
            return args[0]
        res = args[0] & args[1]
        for arg in args[2:]:
            res = res & arg
        return res

    def saturation(self, other):
        """Saturation of ideals (self : other^∞), returns both saturation ideal and saturation index."""
        singular_commands = ["LIB \"elim.lib\";"
                             f"ring r = {self.ring};",
                             f"ideal i = {self};",
                             f"ideal j = {other};",
                             f"list k_and_s = {'sat' if Singular_version <= Version('4.3.1') else 'sat_with_exp'}(i, j);",
                             "print(k_and_s[1]);",
                             "print(\"sat index\"); print(k_and_s[2]);",
                             "$"]
        output = execute_singular_command(singular_commands)
        cls, ring = self.__class__, self.ring
        return (cls(ring, [entry.replace(',', '') for entry in output.split("sat index\n")[0].split("\n") if entry]),
                int(output.split("sat index\n")[1]))

    def __floordiv__(self, other):
        """Saturation of ideals (self : other^∞), returns only the ideal."""
        sat_ideal, _ = self.saturation(other)
        return sat_ideal

    def saturation_index(self, other):
        """Saturation of ideals (self : other^∞), returns only the saturation index."""
        _, sat_index = self.saturation(other)
        return sat_index

    def __truediv__(self, other):
        """Quotient of ideals."""
        singular_commands = [f"ring r = {self.ring};",
                             f"ideal i = {self};",
                             f"ideal j = {other};",
                             "ideal k = quotient(i, j);",
                             "print(k);",
                             "$"]
        output = execute_singular_command(singular_commands)
        output = [line.replace(",", "") for line in output.split("\n")]
        cls, ring = self.__class__, self.ring
        return cls(ring, output)

    def reduce(self, other):
        """Remainder of division, i.e. reduction."""
        if isinstance(other, Ideal):
            singular_commands = [f"ring r = {self.ring};",
                                 f"ideal gb = {','.join(self.groebner_basis)};",
                                 f"ideal i = {other};",
                                 "ideal j = reduce(i, gb);",
                                 "print(j);",
                                 "$"]
            output = execute_singular_command(singular_commands)
            output = output.replace("\n", "").split(",")
            cls, ring = self.__class__, self.ring
            return cls(ring, output)
        else:
            singular_commands = [f"ring r = {self.ring};",
                                 f"ideal gb = {','.join(self.groebner_basis)};",
                                 f"poly f = {other};",
                                 "poly g = reduce(f, gb);",
                                 "print(g);",
                                 "$"]
            output = execute_singular_command(singular_commands)
            return output

    def __str__(self):
        return ",".join(map(str, self.generators))

    def __repr__(self):
        return "Ideal over Ring\n    " + repr(self.ring) + "\ngenerated by\n    " + str(self)

    def to_qring(self, other):
        if self._dim is not None:
            dim_in_full_ring = self._dim
            self.dim_in_full_ring = dim_in_full_ring
        self.delete_cached_properties()
        qring = QuotientRing(self.ring, other)
        self.ring = qring

    def to_full_ring(self):
        if isinstance(self.ring, QuotientRing):
            self.delete_cached_properties()
            self.generators = (self.ring.ideal + self).generators
            self.ring = self.ring.ideal.ring
            if hasattr(self, 'dim_in_full_ring'):
                self.dim = self.dim_in_full_ring
                delattr(self, 'dim_in_full_ring')

    def _get_cached_properties_names(self):
        return [name for name, value in inspect.getmembers(self.__class__)
                if isinstance(value, functools.cached_property) if name in self.__dict__.keys()]

    def delete_cached_properties(self):
        for cached_property in self._get_cached_properties_names():
            delattr(self, cached_property)
        self._dim = None
        self._indepSets = None

    def generators_eval(self, **kwargs):
        return [eval(generator.replace("^", "**"), kwargs) for generator in self.generators]

    @property
    def is_unit_ideal(self):
        return self == Ideal(self.ring, ('1', ))


def reduce(poly, ideal):
    return ideal.reduce(poly)


def monomial_to_exponents(variables, monomial):
    """Converts a monomial in the variables of a polynomial ring into a numpy.array of exponents."""
    exps = numpy.zeros(len(variables), dtype=int)
    variables = list(map(str, variables))  # in case sympy symbols are passsed
    _split_monomial = list(filter(lambda x: x != '', monomial.replace("^", "**").split("*")))
    split_monomial = []
    for entry in _split_monomial:  # need to fix Singular notation for exponents (e.g. if x2 is not a ring variable then x^2 is written as x2)
        if not entry.isdigit() and entry not in variables:
            exp = re.findall(r"(\d+)", entry)[0]
            split_monomial += [entry.replace(exp, ""), exp]
        else:
            split_monomial += [entry]
    for i, ientry in enumerate(split_monomial):
        if not ientry.isdigit():
            if i < len(split_monomial) - 1 and split_monomial[i + 1].isdigit():
                exp = int(split_monomial[i + 1])
            else:
                exp = 1
            exps[variables.index(ientry)] += exp
    return exps
