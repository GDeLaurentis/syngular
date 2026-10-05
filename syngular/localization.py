import re
import sympy

from .ring import Ring
from .ideal import Ideal
from .tools import execute_singular_command


class PrincipalLocalization:
    """Principal localization of a quotient ring at a single element and with a chosen elimination ordering."""

    def __init__(self, qring, element, dependent_variables=(), check_polynomial=False):
        """Construct and cache the reduced Gröbner basis defining the localization."""
        self.qring = qring
        if isinstance(element, (list, tuple)) and all(isinstance(entry, (list, tuple)) and len(entry) == 2 for entry in element):
            self.element_factors = tuple((sympy.sympify(factor), power) for factor, power in element)
        else:
            self.element_factors = ((sympy.sympify(element), 1),)
        self.element = sympy.prod(factor ** power for factor, power in self.element_factors)
        self.dependent_variables = tuple(map(str, dependent_variables))

        if "T" in map(str, qring.ideal.ring.variables):
            raise ValueError("The localization variable T is already present in the ring.")

        remaining_variables = [str(variable) for variable in qring.ideal.ring.variables if str(variable) not in self.dependent_variables]
        self.variables = self.dependent_variables + ("T",) + tuple(remaining_variables)

        self.ring = Ring(qring.ideal.ring.field, self.variables, "lp")
        self.ideal = Ideal(self.ring, tuple(map(str, qring.ideal.generators)) + (f"T*({str(self.element).replace('**', '^')})-1",))
        self.groebner_basis = self.ideal.get_groebner_basis(reduced=True, algorithm="groebner")

        dependent_symbols = set(map(sympy.Symbol, self.dependent_variables))
        excluded_symbols = dependent_symbols | {sympy.Symbol("T")}
        groebner_basis_sympy = [sympy.sympify(entry.replace("^", "**")) for entry in self.groebner_basis]
        self.relations = [entry for entry in groebner_basis_sympy if not entry.free_symbols & excluded_symbols]

        if check_polynomial and self.relations:
            raise RuntimeError(f"Principal localization is not a localization of a polynomial ring. "
                               f"Relations remain among the retained variables:\n{self.relations}")

        if check_polynomial:
            for variable in self.dependent_variables:
                remainder = self.ideal.reduce(variable)
                remainder = sympy.sympify(str(remainder).replace("^", "**"))
                if remainder.free_symbols & dependent_symbols:
                    raise RuntimeError(f"Principal localization did not eliminate {variable}.")

    @property
    def is_polynomial_localization(self):
        return not self.relations

    def __str__(self):
        return f"{self.qring}[1/({self.element})]"

    def __repr__(self):
        return f"PrincipalLocalization({self.qring!r}, {str(self.element)!r}, dependent_variables={self.dependent_variables!r})"

    def __call__(self, polynomial, timeout="default"):
        """Map a polynomial to the localization and return numerator and denominator power."""
        singular_commands = [
            f"ring r = {self.ring};",
            f"poly f = {str(polynomial).replace('**', '^')};",
            f"poly D = {str(self.element).replace('**', '^')};",
            f"ideal G = {','.join(self.groebner_basis)};",
            'attrib(G,"isSB",1);',
            "f = reduce(f,G);",
            "matrix C = coeffs(f,T);",
            "int k = nrows(C)-1;",
            "poly num = 0;",
            "int i;",
            "for (i=0; i<=k; i=i+1)",
            "{",
            "    num = num+C[i+1,1]*D^(k-i);",
            "}",
            'print("__NUM__");',
            "print(num);",
            'print("__POWER__");',
            "print(k);",
            'print("__END__");',
            "$",
        ]

        output = execute_singular_command(singular_commands, timeout=timeout)
        numerator = "".join(output.split("__NUM__", 1)[1].split("__POWER__", 1)[0].splitlines()).strip()
        denominator_power = int(output.split("__POWER__", 1)[1].split("__END__", 1)[0].strip())

        symbols = set(re.findall(r"[A-Za-z_][A-Za-z0-9_]*", numerator))
        dependent_symbols = set(self.dependent_variables)
        if symbols & dependent_symbols:
            raise RuntimeError(f"Principal localization did not eliminate the requested variables: {symbols & dependent_symbols}")

        return numerator, denominator_power
