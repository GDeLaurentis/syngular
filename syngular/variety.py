import functools
import mpmath
import numpy
import random
import re
import sympy
import syngular

from copy import copy, deepcopy
from pyadic import PAdic, ModP

from .tools import RootNotInFieldError, RootPrecisionError

# this fixes a weird bug where sympy does not respect precision even if mpmath.mp.dps precision is set
# (sympy seems to use mpmath as backhand)
equation = sympy.sympify(f"x - 1.{'0' * 290}1")
sympy.nroots(equation, n=300, maxsteps=500)


def retry_to_find_root(max_tries=100):
    def retry_to_find_root_decorator(func):
        @functools.wraps(func)
        def wrapper(self, field, base_point={}, directions=[], valuations=tuple(), indepSetNbr=None, seed=None, verbose=False):
            if base_point != {} and indepSetNbr is not None:
                return func(self, field, base_point=base_point, directions=directions, valuations=valuations,
                            indepSetNbr=indepSetNbr, seed=seed, verbose=verbose)
            else:
                for try_nbr in range(max_tries):
                    try:
                        res = func(self, field, base_point=base_point, directions=directions, valuations=valuations,
                                   indepSetNbr=indepSetNbr, seed=seed, verbose=verbose)
                        break
                    except (RootNotInFieldError, RootPrecisionError) as e:
                        if try_nbr != max_tries - 1:
                            if verbose:
                                print(f"Caught {type(e).__name__}, retrying...")
                            if seed is not None:  # maintain pseudo-randomness, but change seed, else retring has no effect.
                                random.seed(seed)
                                seed += random.randint(10 ** 5, 10**6)
                            continue
                        else:
                            raise type(e)(f"Could not find a root in the field {field} after {max_tries} attempts. {e}")
                return res
        return wrapper
    return retry_to_find_root_decorator


class Variety_of_Ideal:

    @retry_to_find_root(max_tries=100)
    def point_on_variety(self, field, base_point={}, directions=[], valuations=tuple(), indepSetNbr=None, seed=None, verbose=False):
        """Generate a representative point on or close to the variety associated to this ideal.
        The point is 'valuations' away from the exact variety, in the directions specified by 'directions'.
        If 'directions' are not provided, pick the first n=codim simplest generators from 'self'.
        If the ideal is not prime, an irreducible branch will be picked at random."""

        from .ideal import Ideal
        from .qring import QuotientRing

        assert all([valuation > 0 for valuation in valuations])
        self = deepcopy(self)   # don't modify self within this function

        random.seed(seed)

        # handle directions, i.e. the generators of the sub-ideal of maximal codimension
        if directions == []:
            directions = sorted(self.generators, key=lambda x: len(x))[:self.codim]
            assert Ideal(self.ring, directions).codim == self.codim
        # extra directions in qring
        if isinstance(self.ring, QuotientRing):
            directions += self.ring.ideal.generators

        # handle valuations - if valuations == tuple() return a point exactly on the variety
        if field.name == "padic":
            prime, iterations = field.characteristic, field.digits
            if valuations == tuple():
                valuations = tuple(field.digits for poly in directions)
        elif field.name == "finite field":
            prime, iterations = field.characteristic, 1
        elif field.name == "mpc":
            prime, iterations = None, 1 if valuations == tuple() else 2
        # extra valuations in qring
        if isinstance(self.ring, QuotientRing) and field.name == "padic":
            valuations = valuations + (field.digits, ) * len(self.ring.ideal.generators)
        elif isinstance(self.ring, QuotientRing) and field.name == "mpc":
            valuations = valuations + (0, ) * len(self.ring.ideal.generators)

        if verbose:
            print("Directions, valuations:", directions, valuations)

        # if in qring, go to full ring
        if isinstance(self.ring, QuotientRing):
            self.to_full_ring()

        # handle independent set
        indepSets = self.indepSets
        if verbose:
            print("Codimensions:", set(indepSet.count(0) for indepSet in indepSets))
            print("Number of indepSets:", len(indepSets))
        if indepSetNbr is None:
            indepSetNbr = random.randint(0, len(indepSets) - 1)
        indepSet = indepSets[indepSetNbr]
        indepSymbols = tuple([symbol for i, symbol in enumerate(self.ring.variables) if indepSet[i] == 1])
        depSymbols = tuple([symbol for i, symbol in enumerate(self.ring.variables) if indepSet[i] == 0])
        if verbose:
            print("Chosen indepSet:", indepSet, "indepSymbols:", indepSymbols, "depSymbols:", depSymbols)

        # handle the base point, i.e. the values of the independent variables
        if base_point == {}:
            base_point = {indepSymbol: field.random_element() for indepSymbol in indepSymbols}
        else:
            base_point = {sympy.symbols(str(key)): val for key, val in base_point.items()}
        base_point |= {depSymbol: depSymbol for depSymbol in depSymbols}

        oSemiNumericalIdeal = self._semi_numerical_slice(field, directions, valuations, base_point, depSymbols, verbose=False, iteration=0)

        # print(repr(oSemiNumericalIdeal))

        for iteration in range(iterations):

            if verbose:
                print(f"\nAt iteration {iteration}")
                # print(f"base point {base_point}")
                # print(repr(oSemiNumericalIdeal), oSemiNumericalIdeal.primary_decomposition, oSemiNumericalIdeal.groebner_basis, len(oSemiNumericalIdeal.groebner_basis))

            syngular.DEGBOUND.set(0)

            if prime is None and oSemiNumericalIdeal.dim == -1:  # this is the ideal generated by '1'
                raise RootPrecisionError
            assert oSemiNumericalIdeal.dim == 0

            root_dicts = lex_groebner_solve(oSemiNumericalIdeal.groebner_basis, prime=prime)
            check_solutions(oSemiNumericalIdeal.groebner_basis, root_dicts, prime=prime)

            try:
                root_dict = root_dicts[0]
            except IndexError:
                if not field.is_algebraically_closed:
                    raise RootNotInFieldError(f"Got root_dicts: {root_dicts}, for lex Groebner basis:\n{oSemiNumericalIdeal.groebner_basis}.")
                else:
                    raise IndexError(f"Got root_dicts: {root_dicts}, for lex Groebner basis:\n{oSemiNumericalIdeal.groebner_basis}.")
            root_dict = {key: root_dict[key] for key in root_dict.keys() if key not in indepSymbols}

            if iteration < iterations - 1:  # happens only with p-adics
                for key in root_dict.keys():
                    root_dict[key] = root_dict[key] + (prime if prime is not None else 1) * key

            # print("root_dict:, root_dict)

            update_point_dict(base_point, root_dict, field)

            if iteration < iterations - 1:
                if prime is not None:
                    valuations = [valuation - 1 for valuation in valuations]
                    valuations_at_directions = list(filter(lambda x: x[1] > 0, zip(directions, valuations)))
                else:
                    valuations_at_directions = list(zip(directions, valuations))
                if valuations_at_directions != []:
                    directions, valuations = zip(*valuations_at_directions)
                else:
                    directions, valuations = (), ()
                if verbose:
                    print("New directions, valuations:", directions, valuations)

                oSemiNumericalIdeal = self._semi_numerical_slice(field, directions, valuations, base_point, depSymbols, verbose=False, iteration=iteration + 1)

                if len(oSemiNumericalIdeal.indepSets) == 1 and numpy.all(numpy.array(oSemiNumericalIdeal.indepSets) == 0):
                    continue

                currentIndepSet = oSemiNumericalIdeal.indepSets[0]
                newIndepSymbols = tuple([symbol for i, symbol in enumerate(depSymbols) if currentIndepSet[i] == 1])

                if verbose:
                    print("New independent symbols:", newIndepSymbols)

                if newIndepSymbols != tuple():  # this happens only with padics - codimension during iterative lift may drop

                    rand_dict = {newIndepSymbol: random.randrange(1, field.characteristic ** (field.digits - iteration)) for newIndepSymbol in newIndepSymbols}

                    update_point_dict(base_point, rand_dict, field)

                    depSymbols = tuple(symbol for symbol in depSymbols if symbol not in newIndepSymbols)

                    if depSymbols == tuple():  # no more equations to solve, terminate early
                        break

                    oSemiNumericalIdeal = self._semi_numerical_slice(field, directions, valuations, base_point, depSymbols, verbose=False, iteration=iteration + 1)

        if field.name == "padic":
            for key, val in base_point.items():
                if not isinstance(val, PAdic):
                    base_point[key] = PAdic(base_point[key], field.characteristic, field.digits)

        return {str(key): val for key, val in base_point.items()}

    def _semi_numerical_slice(self, field, directions, valuations, base_point, depSymbols, verbose=False, iteration=0):
        """Helper function for point_on_variety. Uses the values in 'base_point' to return a new ideal of lower dimension.
        The potentially perturbed slice of the origial ideal 'self'."""

        from .ideal import Ideal
        from .ring import Ring

        # on subsequent iterations, switch to an ideal of maximal codimension (potentiallty losing branch information), beacuse:
        # (floats) need to append perturbations to equations; (padics) may need to solve less equations; (finite field) iteration > 0 does not make sense, perturbation impossible.
        if iteration > 0:
            generators = list(copy(directions))
            if field.characteristic == 0:
                for i, valuation in enumerate(valuations):
                    generators[i] = str(generators[i]) + "-" + str(sympy.sympify(valuation))
        else:
            generators = copy(self.generators)

        generators = sympy.sympify(generators)
        generators = [sympy.expand(generator.subs(base_point)) for generator in generators]
        generators = list(filter(lambda x: x != 0, generators))

        if field.characteristic == 0:
            generators = [str(generator) for generator in generators]
        else:
            def field_division(match):
                """Patch division as in Fp or Qp, else the polynomial coefficients may not be explicitly in the field."""
                integer_numerator = match.group(1)
                symbolic_monomial = match.group(2)
                if symbolic_monomial is None:
                    symbolic_monomial = ''
                integer_denominator = match.group(3)
                # print(integer_numerator, symbolic_monomial, integer_denominator)
                if field.name == "padic":
                    quotient = PAdic(int(integer_numerator), field.characteristic, field.digits) / int(integer_denominator)
                    return str(int(quotient) * field.characteristic ** quotient.n) + symbolic_monomial
                elif field.name == "finite field":
                    return str(int(ModP(int(integer_numerator), field.characteristic) / int(integer_denominator))) + symbolic_monomial
                else:
                    raise Exception("field not understood while patching sympy division.")
            generators = [re.sub(r"(\d+)(\*[^+-]+)*\/(\d+)", field_division, str(generator)) for generator in generators]
            generators = [re.sub(r"(?<![a-z])(\d+)",
                                 lambda match: str(int(match.group(1)) // field.characteristic ** iteration % field.characteristic),
                                 str(generator)) for generator in generators]
            generators = sympy.sympify(generators)
            generators = list(filter(lambda x: x != 0, generators))

        oZeroDimIdeal = Ideal(Ring(field.singular_notation, depSymbols, "lp"), generators)

        return oZeroDimIdeal


def update_point_dict(base_point_dict, new_vals_dict, field):
    for key in base_point_dict:
        if key in new_vals_dict and base_point_dict[key] == key:
            base_point_dict[key] = new_vals_dict[key]
        elif key in new_vals_dict:
            base_point_dict[key] = base_point_dict[key].subs(new_vals_dict)
        if hasattr(base_point_dict[key], "free_symbols") and base_point_dict[key].free_symbols == set() and field.name == "mpc":
            base_point_dict[key] = mpmath.mpc(base_point_dict[key])
    return base_point_dict


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #
# Univariate solvers
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #


def univariate_floating_point_solver(equation, root_dict):
    """Returns all possible solutions of 'equation' over arbitrary precision complex numbers."""
    equation = sympy.sympify(equation).subs(root_dict)
    free_symbols = list(equation.free_symbols)
    assert len(free_symbols) == 1
    symbol = free_symbols[0]
    solutions = list(map(mpmath.mpc, sympy.nroots(equation, n=300, maxsteps=500)))  # mpmath.polyroots is faster, but the parsing is more complicated
    return update_root_dict(symbol, solutions, root_dict)


def univariate_finite_field_solver(equation, root_dict, prime):
    """Returns all possible solutions of 'equation' over a finite field of cardinality 'prime'.
       If already satisfied returns True, if no solution exists returns False."""
    equation = sympy.sympify(equation).subs(root_dict)
    if isinstance(equation, sympy.core.numbers.Integer) and equation % prime == 0:
        return True
    free_symbols = list(equation.free_symbols)
    if len(free_symbols) < 1:
        return False
    if len(free_symbols) > 1:
        raise Exception("Too many free parameters.")
    symbol = free_symbols[0]
    equation = sympy.poly(equation, modulus=prime)
    if equation == 0:
        return True
    pre_factor, factors = sympy.factor_list(sympy.factor(equation, modulus=prime))
    factors = [factor[0] for factor in factors]
    if pre_factor % prime == 0:
        return True
    linear_factors = [factor for factor in factors if sympy.diff(factor, symbol) == 1]
    if linear_factors == []:
        return False
    solutions = [ModP(int(sympy.solve(factor)[0]), prime) for factor in linear_factors]
    return update_root_dict(symbol, solutions, root_dict)


def update_root_dict(symbol, solutions, root_dict):
    """Given solutions and root_dict returns updated root_dicts."""
    root_dicts = [deepcopy(root_dict)]
    root_dicts[0].update({symbol: solutions[0]})
    for solution in solutions[1:]:
        new_root_dict = deepcopy(root_dict)
        new_root_dict.update({symbol: solution})
        root_dicts.append(new_root_dict)
    return root_dicts


def lex_groebner_solve(equations, prime=None):
    """Returns the variety corresponding to a given zero dimensional ideal in lexicographic groebner basis form.
       The variety take the form of a list of dictionaries for the possible values of the variables."""
    root_dicts = [{}]
    for equation in equations:
        temp_dicts = []
        for i, root_dict in enumerate(root_dicts):
            if prime is None:
                sols = univariate_floating_point_solver(equation, root_dict)
            else:
                sols = univariate_finite_field_solver(equation, root_dict, prime)
            if sols is True:
                temp_dicts += [root_dict]
            elif sols is False:
                continue
            else:
                temp_dicts += sols
        else:
            root_dicts = temp_dicts
    return root_dicts


def check_solutions(equations, root_dicts, prime=None):
    """Checks that all solutions in root_dicts solve the equations."""
    for root_dict in root_dicts:
        check = []
        for equation in equations:
            check += [sympy.sympify(equation)]
            for key, value in root_dict.items():
                if prime is None:
                    check[-1] = sympy.simplify(check[-1].subs(key, value))
                else:
                    check[-1] = check[-1].subs(key, value) % prime
        if prime is None:
            assert all([numpy.isclose(complex(entry), 0) for entry in check])
        else:
            assert all([entry == 0 for entry in check])
