import functools
import re
import sympy
import numpy
import inspect

from .tools import execute_singular_command
from .ring import Ring
from .qring import QuotientRing
from .ideal_algorithms import Ideal_Algorithms


class Ideal(Ideal_Algorithms, object):

    def __init__(self, ring, generators):
        if not isinstance(ring, Ring) or not isinstance(generators, list):
            raise Exception("Invalid Ideal initialisation.")
        self.ring = ring
        self.generators = generators
        self.test_valid_ideal()

    def test_valid_ideal(self):
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
        return hash(", ".join(self.groebner_basis)) + hash(self.ring)

    def __eq__(self, other):
        return self.groebner_basis == other.groebner_basis

    def squoosh(self):
        self.generators = self.minbase

    @functools.cached_property
    def dim(self):
        singular_commands = [f"ring r = {self.ring};",
                             f"ideal gb = {','.join(self.groebner_basis)};",
                             "print(dim(gb));",
                             "$"]
        output = execute_singular_command(singular_commands)
        return int(output)

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

    @functools.cached_property
    def indepSet(self):
        singular_commands = [f"ring r = {self.ring};",
                             f"ideal gb = {','.join(self.groebner_basis)};",
                             "print(indepSet(gb));",
                             "$"]
        output = execute_singular_command(singular_commands)
        return tuple(map(int, output.split(",\n")))

    @functools.cached_property
    def indepSets(self):
        singular_commands = [f"ring r = {self.ring};",
                             # f"ideal gb = {','.join(self.groebner_basis)};",   # this breaks singular variety construction, especially with mpcs
                             f"ideal i = {self};",
                             "ideal gb = groebner(i);",
                             "print(indepSet(gb, 1));",
                             "$"]
        output = execute_singular_command(singular_commands)
        if output == 'empty list':
            return []
        indepSets = [tuple(map(int, line.replace(" ", "").split(","))) for line in output.split("\n") if ":" not in line]
        return indepSets

    @functools.cached_property
    def groebner_basis(self):
        singular_commands = ["option(redSB);",
                             f"ring r = {self.ring};",
                             f"ideal i = {self};",
                             "ideal gb = groebner(i);",
                             "print(gb);",
                             "$"]
        singular_command = "\n".join(singular_commands)
        output = execute_singular_command(singular_command)
        output = [line.replace(",", "") for line in output.split("\n")]
        return output

    @functools.cached_property
    def leadGBmonomials(self):
        """Gives the leading monomials of the Groebner basis polynomials."""
        singular_commands = [f"ring r = {self.ring};",
                             f"ideal gb = {','.join(self.groebner_basis)};",
                             "def lts = lead(gb);",
                             "print(lts);"]
        output = execute_singular_command(singular_commands)
        # print(output)
        output = [line.replace(",", "") for line in output.split("\n")]
        return output

    @functools.cached_property
    def minbase(self):
        singular_commands = [f"ring r = {self.ring};",
                             f"ideal i = {self};",
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
        ring = Ring(ring.field, ring.variables[:var_range.start] + ring.variables[var_range.stop:], ring.ordering)
        return cls(ring, output)

    def __contains__(self, other):
        """Implements ideal membership."""
        if isinstance(other, Ideal):
            return self.__ideal_contains__(other)
        elif isinstance(other, sympy.Expr) or isinstance(other, str):
            return self.__poly_contains__(other)

    def __ideal_contains__(self, other):
        qIdeal = self / other
        if list(map(str, qIdeal.generators)) == ['1']:
            return True
        else:
            return False

    def __poly_contains__(self, other):
        if isinstance(other, sympy.Expr):
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
        self.delete_cached_properties()
        qring = QuotientRing(self.ring, other)
        self.ring = qring

    def to_full_ring(self):
        if isinstance(self.ring, QuotientRing):
            self.delete_cached_properties()
            self.generators = (self.ring.ideal + self).generators
            self.ring = self.ring.ideal.ring

    def _get_cached_properties_names(self):
        return [name for name, value in inspect.getmembers(self.__class__) if isinstance(value, functools.cached_property) if name in self.__dict__.keys()]

    def delete_cached_properties(self):
        for cached_property in self._get_cached_properties_names():
            delattr(self, cached_property)


def reduce(poly, ideal):
    return ideal.reduce(poly)


def monomial_to_exponents(variables, monomial):
    """Converts a monomial in the variables of a polynomial ring into a numpy.array of exponents."""
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
