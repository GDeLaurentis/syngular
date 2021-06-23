import functools
import re
import sympy

from .tools import execute_singular_command
from .ring import Ring
from .qring import QuotientRing


class Ideal(object):

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
    def indepSets(self):
        singular_commands = [f"ring r = {self.ring};",
                             f"ideal i = {self};",
                             "ideal gb = groebner(i);",
                             "print(indepSet(gb, 1));",
                             "$"]
        output = execute_singular_command(singular_commands)
        if output == 'empty list':
            return []
        indepSets = [tuple(map(int, line.replace(" ", "").split(","))) for line in output.split("\n") if":" not in line]
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
    def minbase(self):
        singular_commands = [f"ring r = {self.ring};",
                             f"ideal i = {self};",
                             "print(minbase(i));",
                             "$"]
        output = execute_singular_command(singular_commands)
        output = [line.replace(",", "") for line in output.split("\n")]
        return output

    @functools.cached_property
    def primary_decomposition(self):
        singular_commands = ["LIB \"primdec.lib\";",
                             f"ring r = {self.ring};",
                             f"ideal i = {self};",
                             "ideal gb = groebner(i);",
                             "def pr = primdecGTZ(gb);",
                             "print(pr);"]
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
                             f"ideal i = {self};",
                             "ideal gb = groebner(i);",
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

    def __str__(self):
        return ",".join(map(str, self.generators))

    def __repr__(self):
        return "Ideal over Ring\n    " + repr(self.ring) + "\ngenerated by\n    " + str(self)
