import sympy
import syngular
import warnings

from pycoretools import flatten

from .tools import execute_singular_command


class Ring(object):

    def __init__(self, field, variables, ordering):
        self.field = field
        self.variables = variables
        self.ordering = ordering
        self.test_valid_ring()

    def test_valid_ring(self):
        singular_commands = [f"ring r = {Ring.__str__(self)};",
                             "print(r);"
                             "$"]
        execute_singular_command(singular_commands)

    def __hash__(self):
        return hash(str(self))

    @property
    def field(self):
        return self._field

    @field.setter
    def field(self, field):
        if isinstance(field, tuple) and isinstance(field[0], str):
            field = sympy.symbols(field)
        self._field = field

    @property
    def variables(self):
        return self._variables

    @variables.setter
    def variables(self, variables):
        if not isinstance(variables, tuple):
            raise Exception("Ring variables must be tuple")
        if len(variables) == 0:
            raise Exception("Ring must be defined over at least one variable")
        if isinstance(variables[0], str):
            variables = sympy.symbols(variables)
        self._variables = variables

    @property
    def ordering(self):
        """Monomial ordering"""
        return self._ordering

    @ordering.setter
    def ordering(self, ordering):
        if isinstance(ordering, str):
            assert ordering[:2] in ['lp', 'rp', 'dp', 'Dp', 'ls', 'rs', 'ds', 'Ds', 'wp']
            if ordering[1] == 's':
                print("Warning: chosen ordering is not a well-ordering.")
        elif not isinstance(ordering, tuple):
            raise TypeError(f"Ring ordering has to be either string or tuple, received {ordering} of type {type(ordering)}")
        self._ordering = ordering

    def __str__(self):
        string = f"""{self.field}, {str(self.variables).replace(",)", ")")}, {str(self.ordering).replace("'", "")}"""
        string += ";\nshort=0"
        if syngular.DEGBOUND != 0:
            string += f";\ndegBound = {syngular.DEGBOUND}"
        return string

    def __repr__(self):
        return str(self)

    def __eq__(self, other):
        assert isinstance(self, Ring) and isinstance(other, Ring)
        return self.field == other.field and self.variables == other.variables and self.ordering == other.ordering

    def exact_divide(self, dividend, divisor, timeout="default"):
        """Exactly divide dividend by divisor in the ring, raising ArithmeticError if not divisible."""
        singular_commands = [
            f"ring r = {self};",
            f"poly f = {str(dividend).replace('**', '^')};",
            f"poly g = {str(divisor).replace('**', '^')};",
            "ideal I = g;",
            "ideal F = f;",
            "matrix T = lift(I,F);",
            "ideal Q = std(0);",
            "poly s = reduce(T[1,1],Q);",
            "poly t = reduce(f-g*s,Q);",
            'print("__QUOTIENT__");',
            "print(s);",
            'print("__REMAINDER__");',
            "print(t);",
            'print("__END__");',
            "$",
        ]

        output = execute_singular_command(singular_commands, timeout=timeout)
        quotient = "".join(output.split("__QUOTIENT__", 1)[1].split("__REMAINDER__", 1)[0].splitlines()).strip()
        remainder = "".join(output.split("__REMAINDER__", 1)[1].split("__END__", 1)[0].splitlines()).strip()

        if remainder != "0":
            raise ArithmeticError(f"{dividend} is not exactly divisible by {divisor} in {self}.")

        return quotient

    def factorize(self, polynomial, timeout="default"):
        """Factorize a polynomial in the ambient polynomial ring, returning (factor, multiplicity) pairs."""
        from .qring import QRing

        if isinstance(self, QRing):
            warnings.warn("Factorization in a quotient ring is performed in the ambient polynomial ring,"
                          " ignoring equivalence relations.", stacklevel=2)

        singular_commands = [
            f"ring r = {self};",
            f"poly f = {str(polynomial).replace('**', '^')};",
            "list L = factorize(f);",
            "ideal F = L[1];",
            "intvec E = L[2];",
            "int i;",
            "for (i=1; i<=size(F); i=i+1)",
            "{",
            '    print("__FACTOR_BEGIN__");',
            "    print(E[i]);",
            '    print("__POLY__");',
            "    print(F[i]);",
            '    print("__FACTOR_END__");',
            "}",
            "$",
        ]

        output = execute_singular_command(singular_commands, timeout=timeout)
        factors = []

        for block in output.split("__FACTOR_BEGIN__")[1:]:
            block = block.split("__FACTOR_END__", 1)[0]
            multiplicity, polynomial = block.split("__POLY__", 1)
            factors.append(("".join(polynomial.splitlines()).strip(), int(multiplicity.strip())))

        return factors

    def zero_ideal(self):
        """Returns the zero ideal ⟨0⟩ in the ring."""
        from .ideal import Ideal
        return Ideal(self, [])

    def unit_ideal(self):
        """Returns the unit ideal ⟨1⟩ in the ring."""
        from .ideal import Ideal
        return Ideal(self, ['1'])

    def random_point(self, field, seed=None):
        """Returns a random numerical point in the given field on the zero ideal of the ring."""
        from .point import RingPoint
        j = self.zero_ideal()
        point_val = j.point_on_variety(field=field, seed=seed)
        return RingPoint(self, field, val=point_val)

    def univariate_slice(ring, field, extra_approximate_constraints=(), indepSet=None, seed=None, verbose=False):
        from .polynomial import Polynomial
        from .point import RingPoint
        from .ideal import Ideal
        from .qring import QRing
        if indepSet is None:
            indepSet = (1, ) * len(ring.variables)
        j = Ideal(ring, extra_approximate_constraints)
        p = j.point_on_variety(
            field=field, seed=seed,
            directions=extra_approximate_constraints,
            valuations=(1, ) * len(extra_approximate_constraints)
        )
        point1 = RingPoint(ring, field, val=p)
        # print('approx gens:', [point1(gen) for gen in extra_approximate_constraints])
        xs = tuple(sympy.symbols([f"x{var}" for var in ring.variables]))
        t = sympy.symbols('t')
        point2partial = {str(x): 0 for i, x in enumerate(xs) if indepSet[i] == 0}
        if isinstance(ring, QRing):
            equations = [sympy.poly(sympy.poly(generator).subs({var: var + t * x for var, x in zip(ring.variables, xs)}).subs(point1).subs(point2partial).expand(),
                                    modulus=field.characteristic ** field.digits)
                         for generator in ring.ideal.generators]
            equations = [entry for entry in flatten([sympy.poly(eq, t).all_coeffs() for eq in equations]) if entry != 0]
            if verbose:
                print("Built q-ring equations:")
                print("[" + ",\n ".join(map(str, equations)) + "]")
        else:
            equations = []

        equations_approximate = [sympy.poly(
            sympy.poly(constraint).subs({var: var + t * x for var, x in zip(ring.variables, xs)}).subs(point1).subs(point2partial).expand(),
            modulus=field.characteristic ** field.digits) for constraint in extra_approximate_constraints
        ]
        # print(equations_approximate)
        equations_approximate = [entry for entry in flatten([sympy.poly(eq, t).all_coeffs() for eq in equations_approximate]) if entry != 0]
        # print(equations_approximate)

        ring2 = Ring(field.characteristic, tuple(x for i, x in enumerate(xs) if indepSet[i] == 1), 'dp')
        ideal2 = Ideal(ring2, list(map(str, equations + equations_approximate)))
        # print(ideal2, ideal2.dim)
        point2 = ideal2.point_on_variety(field, directions=(), seed=None if seed is None else seed + 1, verbose=verbose)
        # Should check whether the ideal is the origin
        point2 = point2partial | point2
        if verbose:
            print("Built shift with line coefficients:")
            print(point1)
            print(point2)

        def univariate_slice(t):
            return {key1: point2[key2] * t + point1[key1] for i, (key1, key2) in enumerate(zip(map(str, ring.variables), map(str, xs)))}

        if isinstance(ring, QRing):
            assert all([Polynomial(generator, field).subs(univariate_slice(field.random())) == 0 for generator in ring.ideal.generators])
        return univariate_slice
