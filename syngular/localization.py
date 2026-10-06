import re
import sympy
import functools

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


def find_principal_localization(qring, independent_variables=None, dependent_variables=None, check_polynomial=True):
    """
    Attempt to find a principal localization yielding a polynomial localization for a suitable coordinate projection.

    For A = F[X]/I, split X = Y ⊔ Z, extend I to I^e ⊂ F(Y)[Z], and compute dim_F(Y)(F(Y)[Z]/I^e). If this dimension is one,
    every z ∈ Z is rational in Y; their denominators determine D such that A[D^-1] ≅ F[Y][D^-1].
    """

    variables = tuple(map(str, qring.ideal.ring.variables))

    # Ensure string format for variables
    def variables_as_strings(values):
        if values is None:
            return None
        if isinstance(values, (str, sympy.Symbol)):
            return (str(values),)
        return tuple(map(str, values))

    independent_variables = variables_as_strings(independent_variables)
    dependent_variables = variables_as_strings(dependent_variables)

    # Complete disjoint split if one or both parts were omitted.
    if independent_variables is None and dependent_variables is None:
        indep_set = qring.ideal.indepSet
        if len(indep_set) != len(variables):
            raise RuntimeError("The independent-set mask does not match the quotient-ring variables.")
        independent_variables = tuple(variable for variable, independent in zip(variables, indep_set) if independent)
        dependent_variables = tuple(variable for variable, independent in zip(variables, indep_set) if not independent)
    elif independent_variables is None:
        dependent_set = set(dependent_variables)
        independent_variables = tuple(variable for variable in variables if variable not in dependent_set)
    elif dependent_variables is None:
        independent_set = set(independent_variables)
        dependent_variables = tuple(variable for variable in variables if variable not in independent_set)

    # Sanity checks on provided independent set Y and dependent set Z
    if len(set(independent_variables)) != len(independent_variables) or len(set(dependent_variables)) != len(dependent_variables):
        raise ValueError("Independent and dependent variables must not contain duplicates.")

    unknown_variables = (set(independent_variables) | set(dependent_variables)) - set(variables)
    if unknown_variables:
        raise ValueError(f"Independent/dependent variables contain variables not found in the quotient ring: {sorted(unknown_variables)}")

    if set(independent_variables) & set(dependent_variables):
        raise ValueError("Independent and dependent variables must be disjoint.")

    if set(independent_variables) | set(dependent_variables) != set(variables):
        raise ValueError("Independent and dependent variables must partition the quotient-ring variables.")

    # Trivial case: no dependent variables, so the quotient must already be a polynomial ring.
    if dependent_variables == ():
        if qring.ideal != qring.poly_ring.zero_ideal():
            raise ValueError("No dependent variables were specified, but the quotient ideal is nonzero; check that the independent variables are algebraically independent.")
        return PrincipalLocalization(qring, sympy.Integer(1), (), check_polynomial)

    field = str(qring.poly_ring.field)
    if independent_variables != ():
        independent_string = ",".join(independent_variables)
        field = field[:-1] + f",{independent_string})" if field.startswith("(") and field.endswith(")") else f"({field},{independent_string})"

    extension_ring = Ring(field, dependent_variables, "dp")
    extension_ideal = Ideal(extension_ring, qring.ideal.generators)

    if extension_ideal.dim != 0 or extension_ideal.vdim != 1:
        return None

    normal_forms = extension_ideal.reduce(Ideal(extension_ring, dependent_variables)).generators

    symbols = {variable: sympy.Symbol(variable) for variable in variables}
    dependent_symbols = {symbols[variable] for variable in dependent_variables}

    denominators = []
    for normal_form in normal_forms:
        normal_form = sympy.cancel(sympy.sympify(str(normal_form).replace("^", "**"), locals=symbols))

        if normal_form.free_symbols & dependent_symbols:
            raise RuntimeError(f"Dependent variables remain after reduction: {normal_form}")

        denominators.append(normal_form.as_numer_denom()[1])

    element = sympy.expand(functools.reduce(sympy.lcm, denominators, sympy.Integer(1)))

    return PrincipalLocalization(qring, element, dependent_variables, check_polynomial)
