import sympy
import re

from random import randint
from copy import deepcopy
from packaging.version import Version

from .tools import execute_singular_command, Singular_version
from .ring import Ring
from .field import Field
from .polynomial import Polynomial
from .settings import TemporarySetting


class Ideal_Algorithms:

    singular_commands_EXT = r"""
        // option(redSB);
        ring r = {r};
        ideal i = {ideal};
        // Block Order X \ U > U in F[X]
        ring r1 = {r1};
        ideal j = imap(r, i);
        ideal gr = groebner(j);
        // to F(U)[X \ U]
        ring r2 = {r2};
        ideal G = imap(r1, gr);
        print(G);
    """

    def extension(self, U, ordering="lp"):
        r"""Extension to K(U)[X \ U]"""
        r = self.ring
        X = r.variables
        XqU = tuple(entry for entry in X if entry not in U)
        if ordering == "lp":
            r1 = Ring('0', XqU + U, 'lp')
            r2 = Ring((sympy.symbols('0'), ) + U, XqU, 'lp')
        elif ordering == "dp":
            r1 = Ring('0', XqU + U, (f'dp({len(XqU)})', f'dp({len(U)})'))
            r2 = Ring((sympy.symbols('0'), ) + U, XqU, 'dp')
        else:
            raise ValueError
        string = execute_singular_command(self.singular_commands_EXT.format(**{"ideal": self, "r": r, "r1": r1, "r2": r2}))
        generators = string.split(",\n")
        Ideal = self.__class__
        I_extended = Ideal(r2, generators)
        return I_extended

    singular_commands_EXTCONT1 = r"""
        LIB "polylib.lib";
        // K(U)[X \ U]
        ring r2 = {r2};
        ideal G = {extended_ideal};
        // List of leadcoefficients
        list l;
        for (int i=1; i<=size(G); i=i+1)
        {{
            l = insert(l, leadcoef(G[i]));
        }}
        ring r = {r};
        list l = imap(r2, l);
        // f-polynomial
        list factors;
        for (i=1; i<=size(l); i=i+1)
        {{
            factors = insert(factors, factorize(l[i], 1));
        }}
        print("factors:"); print(factors);
        $
    """

    def extension_contraction_fpoly(self, U, ordering="lp"):
        """Returns a minimal list of factors for the f-polynomial in CONT and EXTCONT"""
        r = self.ring
        X = r.variables
        XqU = tuple(entry for entry in X if entry not in U)
        if ordering == "lp":
            r1 = Ring('0', XqU + U, 'lp')
            r2 = Ring((sympy.symbols('0'), ) + U, XqU, 'lp')
        elif ordering == "dp":
            r1 = Ring('0', XqU + U, (f'dp({len(XqU)})', f'dp({len(U)})'))
            r2 = Ring((sympy.symbols('0'), ) + U, XqU, 'dp')
        else:
            raise ValueError
        string = execute_singular_command(self.singular_commands_EXTCONT1.format(**{"extended_ideal": self.extension(U, ordering), "r": r, "r1": r1, "r2": r2}))
        f_polys_factors = list(set([row.split("=")[1].replace(" ", "") for row in string.split("\n") if "factors" not in row and "=" in row]))
        f_polys_factors = [entry for entry in f_polys_factors if entry != '1']
        return ['1'] + f_polys_factors

    singular_commands_EXTCONT2 = """
        LIB "polylib.lib";
        // option(redSB);
        ring r2 = {r2};
        ideal G = {extended_ideal};
        ring r = {r};
        ideal Iec = imap(r2, G);
        list factors = {f_polys_factors};
        list Iec_and_s;
        int s = 0;
        for (int i_=1; i_<=size(factors); i_=i_+1)
        {{
            Iec_and_s = {sat}(Iec, ideal(factors[i_]));
            Iec = Iec_and_s[1];
            s = max(s, Iec_and_s[2]);
        }}
        print("saturation index:"); print(s);
        print("extended-contracted ideal:"); print(minbase(Iec));
    """

    def extension_contraction(self, U, ordering="lp"):
        """Computes extension-contraction of self to localization define by indepSet U. Returns a tuple: (saturation index, extended-contracted ideal)."""
        r = self.ring
        X = r.variables
        XqU = tuple(entry for entry in X if entry not in U)
        if ordering == "lp":
            # r1 = Ring('0', XqU + U, 'lp')
            r2 = Ring((sympy.symbols('0'), ) + U, XqU, 'lp')
        elif ordering == "dp":
            # r1 = Ring('0', XqU + U, (f'dp({len(XqU)})', f'dp({len(U)})'))
            r2 = Ring((sympy.symbols('0'), ) + U, XqU, 'dp')
        else:
            raise ValueError
        string = execute_singular_command(self.singular_commands_EXTCONT2.format(**{
            "extended_ideal": self.extension(U, ordering), "f_polys_factors": ','.join(self.extension_contraction_fpoly(U, ordering)),
            "r": r, "r2": r2, "sat": "sat" if Singular_version <= Version('4.3.1') else 'sat_with_exp'}))
        Ideal = self.__class__
        return int(string.split("\n")[1]), Ideal(r, [entry.replace(",", "") for entry in string.split("\n")[3:]])

    def primeTestDLP(self, verbose=False, timeout_fpoly=10, timeout_dim=600,
                     seminumerical_dim_computation=False, nbr_points=100,
                     iterated_degbound_computation=False, projection_number=None,
                     fpoly_number=None, astuple=False):
        """
        Assumes equidimensionality of input ideal.
        Returns True if the ideal is prime, False if it is not.
        If astuple is set to True, return (is_primary, is_prime) bool pair.
        Raises Inconclusive if it can't decide.
        Experimental new feature with iterated_degbound_computation=True,
        may help when ideal is primary and deg-unbounded computation fails.
        """
        # todo: rename projection_number to reflect that it can accept directly an indepSet tuple
        import syngular
        if seminumerical_dim_computation and iterated_degbound_computation:
            raise Exception("Semi-numerical test and iterated degbound test are exclusive.")
        # algorithm works over rings, if in a qring convert to the full ring.
        self = deepcopy(self)
        self.to_full_ring()
        # first step: find zero dimensional projections which are linear in the dependent variables, i.e. are a single point.
        if seminumerical_dim_computation:
            try:
                self.indepSets
                self.dim
            except TimeoutError:
                if verbose:
                    print("indepSet and dim computation timedout - will learn semi-numerically.")
            field = Field("finite field", 2 ** 31 - 1, 1)
            points = []
            for i in range(nbr_points):
                if verbose:
                    print(f"\rGenerating point #{i}        ", end="")
                points += [self.point_on_variety(field, indepSet='force guess', seed=i)]
        if projection_number is not None:
            if isinstance(projection_number, int):
                self.indepSets = self.indepSets[projection_number:projection_number + 1]
            elif isinstance(projection_number, tuple):
                self.indepSets = [projection_number]
        lowest_degree_projection_indepSets = []
        lowest_degree = 9999999
        for i, indepSet in enumerate(self.indepSets):
            if verbose:
                print(f"\r@{i}/{len(self.indepSets)}", end="")
            X = self.ring.variables                                            # all variables: X
            U = tuple(var for is_indep, var in zip(indepSet, X) if is_indep)   # independent variables: U
            XqU = tuple(entry for entry in X if entry not in U)                # dependent variables X / U
            prime = 536870909  # largest prime below 2 ** 29
            projection = {entry: randint(1, prime) for entry in U}             # random values for U
            zeroDimSelf = deepcopy(self)
            zeroDimSelf.generators = [str(sympy.poly(entry.subs(projection), modulus=prime).as_expr())
                                      for entry in sympy.sympify(zeroDimSelf.generators)]
            zeroDimSelf.ring = Ring(str(prime), XqU, 'lp')
            zeroDimSelf.delete_cached_properties()
            degrees = list(map(int, re.findall(r"\^(\d)", "".join(zeroDimSelf.groebner_basis))))
            if degrees == []:
                max_degree = 1
            else:
                max_degree = max(degrees)
            if max_degree < lowest_degree:
                lowest_degree = max_degree
                lowest_degree_projection_indepSets = []
            if max_degree == lowest_degree:
                lowest_degree_projection_indepSets += [indepSet]
            if verbose:
                print(f"\r@{i}/{len(self.indepSets)} projections found: {len(lowest_degree_projection_indepSets)} of degree {lowest_degree}     ", end="")
        if verbose:
            print(f"\rprojections found: {len(lowest_degree_projection_indepSets)} of degree {lowest_degree}             ", end="                    ")
        # for each single-point zero-dimensional projection, get the factors of the f-polynomial in EXTCONT (Becker et al. Proposition 8.96)
        f_polys_factors = []
        for i, indepSet in enumerate(lowest_degree_projection_indepSets):
            if verbose:
                number_of_timedout_fpolys = ["TIMEDOUT" in entry[-1] for entry in f_polys_factors].count(True)
                print(f"\rgathering f-poly factors: @ {i}/{len(lowest_degree_projection_indepSets)} of which {number_of_timedout_fpolys} timedout",
                      end="                                   ")
            # 0 dim slice
            X = self.ring.variables
            U = tuple(var for is_indep, var in zip(indepSet, X) if is_indep)
            try:
                with TemporarySetting("syngular", "TIMEOUT", timeout_fpoly):
                    f_polys_factors += [self.extension_contraction_fpoly(U, "lp")]
            except TimeoutError:
                f_polys_factors += [['- TIMEDOUT - probably very very long ' * 8000]]
                continue
            except Exception as e:
                if verbose:
                    print(f"An exception occurred: {e}")
                f_polys_factors += [['- TIMEDOUT - probably very very long ' * 8000]]
                continue
        # just keep the one with the smallest greatest factor
        max_lengths = [max(map(len, f_poly_factors)) for f_poly_factors in f_polys_factors]
        easiest_projection_indepSet = lowest_degree_projection_indepSets[max_lengths.index(min(max_lengths))]
        easiest_projection_indepSet_index = self.indepSets.index(easiest_projection_indepSet)
        if verbose:
            print(f"\neasiest projection is {easiest_projection_indepSet_index}: {easiest_projection_indepSet} of degree {lowest_degree}")
        smallest_fpoly_factors = f_polys_factors[max_lengths.index(min(max_lengths))]
        if lowest_degree > 1:
            tries_for_irreducibility = 100
            for i in range(tries_for_irreducibility):
                # Warning repeated code from above, use _semi_numerical_slice?
                indepSet = easiest_projection_indepSet
                X = self.ring.variables                                            # all variables: X
                U = tuple(var for is_indep, var in zip(indepSet, X) if is_indep)   # independent variables: U
                XqU = tuple(entry for entry in X if entry not in U)                # dependent variables X / U
                prime = 536870909  # largest prime below 2 ** 29
                projection = {entry: randint(1, prime) for entry in U}             # random values for U
                zeroDimSelf = deepcopy(self)
                zeroDimSelf.generators = [str(sympy.poly(entry.subs(projection), modulus=prime).as_expr())
                                          for entry in sympy.sympify(zeroDimSelf.generators)]
                zeroDimSelf.ring = Ring(str(prime), XqU, 'lp')
                zeroDimSelf.delete_cached_properties()
                zeroDim_prim_dec = zeroDimSelf.primary_decomposition
                if len(zeroDim_prim_dec) == 1:
                    Q, P = zeroDim_prim_dec[0]
                    radical = True if Q == P else False
                    print("\rWas irreducible at the chosen point. Continuing test.", end="                    ")
                    break  # conclusive
                elif verbose:
                    print(f"\rWas reducible at the chosen points, checking again. At try {i}.", end="")
            else:
                raise Inconclusive(f"Likely reducible, since it was reducible on {tries_for_irreducibility} zero dim extensions.")
                # return False if not astuple else (False, False)  # statistical statement, may return false negatives
        else:
            radical = True
        if verbose:
            print(f"The ideal is radical: {radical}")
        if smallest_fpoly_factors == ['- TIMEDOUT - probably very very long ' * 8000]:
            raise Inconclusive("Timedout on fpoly factors gathering.")
        smallest_fpoly_factors = sorted(smallest_fpoly_factors, key=lambda x: len(x))
        if verbose:
            print(f"\r smallest f poly factors ({len(smallest_fpoly_factors)}): complexity {sum(map(len, smallest_fpoly_factors))}",
                  end="                              \n")
        # check that the dimensionality drops when adding each of these factors separately (and hence drops for <ideal, f^s>)
        with TemporarySetting("syngular", "TIMEOUT", timeout_dim):
            self.codim  # just cache it - if seminumerical_dim_computation then it was learnt numerically
            if verbose:
                print(f"Original ideal has codim {self.codim}")
            for i, factor in enumerate(smallest_fpoly_factors):
                if fpoly_number is not None and i < fpoly_number:
                    continue
                if verbose:
                    print(f"\r at factor {i}: {factor}.", end="                                       \n")
                if seminumerical_dim_computation:
                    valuations = [Polynomial(factor, field).subs(point) for point in points]
                    zero_count = sum(1 for valuation in valuations if valuation == 0)
                    if verbose and nbr_points != 0:
                        print(f"{zero_count} zeros out of {len(valuations)} points.")
                    if any([valuation == 0 for valuation in valuations]):
                        return False if not astuple else (False, False)
                    else:
                        pass
                        # continue  # statistical, may return false positives
                X = deepcopy(self)
                X.generators += [factor]
                X.delete_cached_properties()
                # Experimental - Assumes codim w/ deg bound <= true codim.
                # Helps termiante the prime test early, IF the result is True.
                if factor != '1' and seminumerical_dim_computation:
                    X.codim_upper_bound = self.codim + 1
                    X.point_on_variety(field, verbose=False)  # learn the dimension of X numerically - possible improvement: verbosity level
                    assert hasattr(X, '_dim') and isinstance(X._dim, int)
                    if verbose:
                        print(f"Learnt codim of Ideal with f-poly factor, {X.codim}")
                for deg_bound in syngular.DEGBOUNDs * iterated_degbound_computation:
                    with TemporarySetting("syngular", "DEGBOUND", deg_bound):
                        X.delete_cached_properties()
                        if self.codim < X.codim:
                            if verbose:
                                print(f"deg_bound {deg_bound} computation was conclusive.", end="\n")
                            break
                        else:
                            if verbose:
                                print(f"deg_bound {deg_bound} computation was inconclusive.", end="\n")

                else:  # loop completed without encountering a break
                    if verbose and iterated_degbound_computation:
                        print("deg_bound reset to zero. Performing full computation.", end="\n")
                    with TemporarySetting("syngular", "DEGBOUND", 0):
                        if not seminumerical_dim_computation:
                            X.delete_cached_properties()
                        # if self.indepSet.count(0) >= X.indepSet.count(0):   # deprecated, equivalent to next line.
                        if self.codim >= X.codim:
                            return False if not astuple else (False, False)
        return radical if not astuple else (True, radical)

    def test_primality(self, *args, **kwargs):
        """Tests if an ideal is prime. Not to be confused with primary."""
        # eventually check assumptions and decide which test to use
        return self.primeTestDLP(*args, **kwargs)


class Inconclusive(Exception):
    pass
