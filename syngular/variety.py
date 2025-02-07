import functools
import mpmath
import numpy
import random
import re
import sympy
import syngular
import warnings

from copy import copy, deepcopy
from pyadic import ModP

from mpmath.libmp.libhyper import NoConvergence

from .tools import RootNotInFieldError, RootPrecisionError
from .field import Field
from .polynomial import Monomial, Polynomial
from .settings import TemporarySetting

# this fixes a weird bug where sympy does not respect precision even if mpmath.mp.dps precision is set
# (sympy seems to use mpmath as backhand)
equation = sympy.sympify(f"x - 1.{'0' * 290}1")
sympy.nroots(equation, n=300, maxsteps=500)


def retry_to_find_root(max_tries=100):
    def retry_to_find_root_decorator(func):
        @functools.wraps(func)
        def wrapper(self, field, base_point={}, directions=None, valuations=tuple(), indepSetNbr=None, indepSet='guess',
                    seed=None, verbose=False):

            if indepSetNbr is not None and indepSet == 'guess':
                indepSet = indepSetNbr
                warnings.warn(
                    "Key word argument 'indepSetNbr' is deprecated and will be removed in a future version. "
                    "Please use 'indepSet' instead.", DeprecationWarning, stacklevel=2
                )

            if indepSet == 'guess':  # check if indep set can be computed quickly
                if verbose:
                    print("At first try, trying to compute the indepSets. ", end="")
                with TemporarySetting("syngular", "TIMEOUT", 3):
                    try:
                        self.indepSets
                        indepSet = None  # if this suceeds there is no reason to guess
                        print("IndepSet computation is cheap. Will compute it.")
                    except TimeoutError:
                        if verbose:
                            print("Gave up after 3 seconds. Will proceed with guesses.")

            if base_point != {} and indepSet not in [None, 'guess']:
                return func(self, field, base_point=base_point, directions=directions, valuations=valuations,
                            indepSet=indepSet, seed=seed, verbose=verbose)
            else:
                for try_nbr in range(max_tries):
                    try:
                        res = func(self, field, base_point=base_point, directions=directions, valuations=valuations,
                                   indepSet=indepSet, seed=seed, verbose=verbose)
                        break
                    except (RootNotInFieldError, RootPrecisionError, NoConvergence, AssertionError) as e:
                        if try_nbr != max_tries - 1:
                            if verbose:
                                print(f"Caught {type(e).__name__} at try number {try_nbr}, retrying...")
                            if seed is not None:  # maintain pseudo-randomness, but change seed, else retring has no effect.
                                random.seed(seed)
                                seed += random.randint(10 ** 5, 10**6)
                            continue
                        else:
                            raise type(e)(f"Could not find a solution in {field} after {max_tries} attempts. {e}")
                return res
        return wrapper
    return retry_to_find_root_decorator


class Variety_of_Ideal:

    @retry_to_find_root(max_tries=100)
    def point_on_variety(self, field, base_point={}, directions=None, valuations=tuple(), indepSet='guess',
                         seed=None, verbose=False, directions_analytic_check=False):
        """Generate a representative point on or close to the variety associated to this ideal.
        The point is 'valuations' away from the exact variety, in the directions specified by 'directions'.
        If 'directions' are not provided, pick the first n=codim simplest generators from 'self'.
        If the ideal is not prime, an irreducible branch will be picked at random."""

        from .ideal import Ideal
        from .qring import QuotientRing

        assert all([isinstance(valuation, str) or valuation >= 0 for valuation in valuations])
        original_self = self
        self = deepcopy(self)   # don't modify self within this fuanction

        random.seed(seed)

        # do not modify directions, in case re-try is triggered, better the input is identical
        directions = deepcopy(directions)
        # handle directions, i.e. the generators of the sub-ideal of maximal codimension
        if directions is not None and directions != [] and field.name == "finite field":
            if verbose:
                print("Warning: directions passed with a finite field. Discarding them.")
            directions = []
        if directions is None or directions == []:
            directions = []
            if field.name != "finite field":
                if verbose:
                    print("Directions not provided, obtaining them from ideal generators.")
                for poly in sorted(self.generators, key=lambda x: len(x)):
                    if Ideal(self.ring, directions + [poly, ]).codim > Ideal(self.ring, directions).codim:
                        directions += [poly, ]
                    if len(directions) == self.codim:
                        break
                assert Ideal(self.ring, directions).codim == self.codim
        elif directions_analytic_check:
            # make sure the provided directions make sense, given the ideal (they need to belong to it)
            # if this is not triggered, then the check is perfomed numerically later. The analytic check can be expensive.
            for direction in directions:
                if Polynomial(direction, Field("rational", 0, 0)) not in self:
                    raise Exception(f"Invalid direction, {direction} was not in {self}.")
        # make sure provided directions are expanded strings
        for i, direction in enumerate(directions):
            if isinstance(direction, sympy.core.Basic):
                directions[i] = str(sympy.expand(direction))
        # extra directions in qring
        if isinstance(self.ring, QuotientRing) and field.name not in ['finite field', 'Fp']:
            directions += self.ring.ideal.generators

        # handle valuations - if valuations == tuple() return a point exactly on the variety
        if field.name == "padic":
            prime, iterations = field.characteristic, field.digits
            if valuations == tuple():
                valuations = tuple(field.digits for _ in directions)
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
            print("Directions:", directions, )
            print("Valuations:", valuations)

        # if in qring, go to full ring
        if isinstance(self.ring, QuotientRing):
            self.to_full_ring()

        # handle independent set
        if indepSet == 'guess':
            chose_indepSet = self.guess_indep_set()
        elif indepSet is None or isinstance(indepSet, int):
            indepSets = self.indepSets
            if verbose:
                print("Codimensions:", set(indepSet.count(0) for indepSet in indepSets))
                print("Number of indepSets:", len(indepSets))
            chose_indepSet = indepSets[random.randint(0, len(indepSets) - 1) if indepSet is None else indepSet]
        indepSymbols = tuple([str(symbol) for i, symbol in enumerate(self.ring.variables) if chose_indepSet[i] == 1])
        depSymbols = tuple([str(symbol)for i, symbol in enumerate(self.ring.variables) if chose_indepSet[i] == 0])
        if verbose:
            print("Chosen indepSet:", chose_indepSet)
            print("indepSymbols:", indepSymbols)
            print("depSymbols:", depSymbols)
            print(f"{'Guessing codim' if indepSet == 'guess' and self._dim is None else 'Codim'} {len(depSymbols)} variety in {len(chose_indepSet)}-dim space")

        # handle the base point, i.e. the values of the independent variables
        if base_point == {}:
            base_point = {indepSymbol: field.random() for indepSymbol in indepSymbols}
        else:
            base_point = {str(key): field(val) for key, val in base_point.items()}
        base_point |= {depSymbol: Polynomial(depSymbol, field) for depSymbol in depSymbols}

        oSemiNumericalIdeal = self._semi_numerical_slice(field, directions, valuations, base_point, depSymbols, verbose=verbose, iteration=0)

        # print(repr(oSemiNumericalIdeal))

        for iteration in range(iterations):

            if verbose:
                print(f"\nAt iteration {iteration}")
                # print(f"base point {base_point}")
                # print(repr(oSemiNumericalIdeal), oSemiNumericalIdeal.primary_decomposition, oSemiNumericalIdeal.groebner_basis, len(oSemiNumericalIdeal.groebner_basis))

            syngular.DEGBOUND = 0

            # and oSemiNumericalIdeal.dim == -1:  # this is the ideal generated by '1'
            # check it explicitly, unless the .dim property is reverted to use grobner_basis instead of std
            if prime is None and oSemiNumericalIdeal.groebner_basis == ['1']:
                raise RootPrecisionError
            if not oSemiNumericalIdeal.dim == 0:
                if oSemiNumericalIdeal.dim > 0 and indepSet == 'guess':
                    # determines dimension of original ideal in the full ring from that of the semi-numerical slice
                    learnt_dimension = len(self.ring.variables) - len(self.generators) + oSemiNumericalIdeal.dim
                    if isinstance(original_self.ring, QuotientRing):
                        original_self.dim_in_full_ring = learnt_dimension
                    else:
                        original_self.dim = learnt_dimension
                    if verbose:
                        print(f"Determined the actual dimension to be {self.dim}")
                raise AssertionError(f"The dimension of the semi-numerical ideal was {oSemiNumericalIdeal.dim} instead of zero.")

            root_dicts = lex_groebner_solve(oSemiNumericalIdeal.groebner_basis, prime=prime)
            check_solutions(oSemiNumericalIdeal.groebner_basis, root_dicts, field)  # they may be stricter then wanted for mpc.

            try:
                root_dict = root_dicts[0]
            except IndexError:
                if not field.is_algebraically_closed:
                    raise RootNotInFieldError(f"Got root_dicts: {root_dicts}, for lex Groebner basis:\n{oSemiNumericalIdeal.groebner_basis}.")
                else:
                    raise IndexError(f"Got root_dicts: {root_dicts}, for lex Groebner basis:\n{oSemiNumericalIdeal.groebner_basis}.")
            # root_dict = {key: root_dict[key] for key in root_dict.keys() if key not in indepSymbols}

            if iteration < iterations - 1:
                for key in root_dict.keys():
                    if field.name == "padic":
                        root_as_poly = Polynomial(root_dict[key], Field("finite field", field.characteristic, 1))
                        root_as_poly.field = field
                    else:
                        root_as_poly = Polynomial(root_dict[key], field)
                    root_dict[key] = root_as_poly + Polynomial(f"{(prime if prime is not None else 1)} * {key}", field)

            update_point_dict(base_point, root_dict, field)
            # print("updated point:", base_point)

            if iteration == 0 and not directions_analytic_check:
                # instead of analytically checking the directions are consistent with the ideal
                # we check they vanish numerically at the constructed point
                for direction in directions:
                    num_poly = Polynomial(direction, field).subs(base_point).subs({key: 0 for key in depSymbols}).coeffs[0]
                    # print("val:", float(abs(num_poly)))
                    if abs(num_poly) > abs(field.ε):
                        raise Exception(f"Invalid direction, {direction} was not in {self}. Numerical membership check failed.")

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

                    if field.name != "padic":
                        raise Exception(f"Codimension changed while not using padics, the field was {field}. Are you sure number of valuations is correct?")

                    rand_dict = {newIndepSymbol: random.randrange(1, field.characteristic ** (field.digits - iteration)) for newIndepSymbol in newIndepSymbols}

                    update_point_dict(base_point, rand_dict, field)

                    depSymbols = tuple(symbol for symbol in depSymbols if symbol not in newIndepSymbols)

                    if depSymbols == tuple():  # no more equations to solve, terminate early
                        break

                    oSemiNumericalIdeal = self._semi_numerical_slice(field, directions, valuations, base_point, depSymbols, verbose=False, iteration=iteration + 1)

        for key, val in base_point.items():
            if val not in field:
                if not isinstance(base_point[key], Polynomial):
                    raise ValueError(f"{base_point[key]} is {type(base_point[key])} instead of a Polynomial")
                assert len(base_point[key]) == 1 and base_point[key].coeffs_and_monomials[0][1] == Monomial("")
                base_point[key] = base_point[key].coeffs_and_monomials[0][0]

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
                    generators[i] = str(generators[i]) + " + " + str(-field(valuation))
        else:
            generators = copy(self.generators)

        generators = [Polynomial(generator, field=field) for generator in generators]
        generators = [generator.subs(base_point, field=field) for generator in generators]
        generators = list(filter(lambda x: x != 0, generators))

        if field.characteristic == 0:
            generators = [re.sub(r"(\d)j", r"\1*I", str(generator), ) for generator in generators]
        elif field.name == "padic":
            for i, generator in enumerate(generators):
                generators[i] = generator / field.characteristic ** iteration
                generators[i].coeffs = [coeff.as_tuple_from_zero[0] for coeff in generators[i].coeffs]
                generators[i]._field = Field("finite field", field.characteristic, 1)
            generators = list(filter(lambda x: x != 0, generators))

        oZeroDimIdeal = Ideal(Ring(field.singular_notation, depSymbols, "lp"), generators)

        if verbose:
            print("Constructed semi-numerical ideal slice:")
            print(repr(oZeroDimIdeal))

        return oZeroDimIdeal


def update_point_dict(base_point_dict, new_vals_dict, field):
    for key in base_point_dict:
        if key in new_vals_dict and base_point_dict[key] == key:
            base_point_dict[key] = new_vals_dict[key]
        elif key in new_vals_dict:
            base_point_dict[key] = base_point_dict[key].subs(new_vals_dict)
        elif key in list(map(str, new_vals_dict.keys())):
            base_point_dict[key] = base_point_dict[key].subs({str(key): val for key, val in new_vals_dict.items()})
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


def check_solutions(equations, root_dicts, field):
    """Checks that all solutions in root_dicts solve the equations."""
    field = field if field.name not in ["padic", "Qp"] else Field("finite field", field.characteristic, 1)
    for root_dict in root_dicts:
        check = [Polynomial(equation, field).subs(root_dict) for equation in equations]
        assert all([len(entry) == 1 for entry in check])  # check these are elements of the field
        if not all([abs(entry.coeffs[0]) <= field.tollerance for entry in check]):
            if field.characteristic == 0:
                raise RootPrecisionError
            else:
                raise AssertionError
