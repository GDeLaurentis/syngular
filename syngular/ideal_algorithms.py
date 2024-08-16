import sympy
import re

from random import randint
from copy import deepcopy

from .tools import execute_singular_command
from .ring import Ring


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
            Iec_and_s = sat(Iec, ideal(factors[i_]));
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
        string = execute_singular_command(self.singular_commands_EXTCONT2.format(**{"extended_ideal": self.extension(U, ordering),
                                                                                    "f_polys_factors": ','.join(self.extension_contraction_fpoly(U, ordering)), "r": r, "r2": r2}))
        Ideal = self.__class__
        return int(string.split("\n")[1]), Ideal(r, [entry.replace(",", "") for entry in string.split("\n")[3:]])

    def primeTestDLP(self, verbose=False, timeout_fpoly=10, timeout_dim=600,
                     iterated_degbound_computation=False, projection_number=None):
        """Returns True if the ideal is prime, False if it is not. Raises Inconclusive if it can't decide. Assumes equidimensionality of input ideal.
        Experimental new feature with iterated_degbound_computation=True, may help when ideal is prime and deg-unbounded computation fails.
        """
        import syngular
        # algorithm works over rings, if in a qring convert to the full ring.
        self.to_full_ring()
        # first step: find zero dimensional projections which are linear in the dependent variables, i.e. are a single point.
        if projection_number is None:
            linear_projection_indepSetIndices = []
            for i in range(len(self.indepSets)):
                if verbose:
                    print(f"\r @ {i}/{len(self.indepSets)}", end="")
                indepSet = self.indepSets[i]
                X = self.ring.variables                                            # all variables: X
                U = tuple(var for is_indep, var in zip(indepSet, X) if is_indep)   # independent variables: U
                XqU = tuple(entry for entry in X if entry not in U)                # dependent variables X / U
                projection = {entry: randint(1, 2 ** 31 - 19) for entry in U}      # random values for U
                zeroDimSelf = deepcopy(self)
                zeroDimSelf.generators = [str(sympy.poly(entry.subs(projection), modulus=2 ** 31 - 19).as_expr()) for entry in sympy.sympify(zeroDimSelf.generators)]
                zeroDimSelf.ring = Ring(str(2 ** 31 - 19), XqU, 'lp')
                zeroDimSelf.delete_cached_properties()
                degrees = list(map(int, re.findall(r"\^(\d)", "".join(zeroDimSelf.groebner_basis))))
                if degrees == []:
                    max_degree = 1
                else:
                    max_degree = max(degrees)
                if max_degree == 1:
                    linear_projection_indepSetIndices += [i]
                if verbose:
                    print(f"\r @ {i}/{len(self.indepSets)} linear projections found: {len(linear_projection_indepSetIndices)}", end="")
            if verbose:
                print(f"\r linear projections found: {len(linear_projection_indepSetIndices)}", end="                    ")
            if len(linear_projection_indepSetIndices) == 0:
                raise Inconclusive
            # for each single-point zero-dimensional projection, get the factors of the f-polynomial in EXTCONT (Becker et al. Proposition 8.96)
            f_polys_factors = []
            for i, linear_projection_indepSetIndex in enumerate(linear_projection_indepSetIndices):
                if verbose:
                    print(f"\r gathering f-poly factors: @ {i}/{len(linear_projection_indepSetIndices)}", end="                                   ")
                indepSet = self.indepSets[linear_projection_indepSetIndex]
                # 0 dim slice
                X = self.ring.variables
                U = tuple(var for is_indep, var in zip(indepSet, X) if is_indep)
                try:
                    f_polys_factors += [self.extension_contraction_fpoly(U, "lp")]
                except TimeoutError:
                    f_polys_factors += [['- TIMEDOUT - probably very very long ' * 80]]
                    continue
                except Exception as e:
                    if verbose:
                        print(f"An exception occurred: {e}")
                    f_polys_factors += [['- TIMEDOUT - probably very very long ' * 80]]
                    continue
            # just keep the one with the smallest greatest factor
            syngular.TIMEOUT = timeout_fpoly
            max_lengths = [max(map(len, f_poly_factors)) for f_poly_factors in f_polys_factors]
            if verbose:
                print("\n easiest projection:", linear_projection_indepSetIndices[max_lengths.index(min(max_lengths))])
            smallest_fpoly_factors = f_polys_factors[max_lengths.index(min(max_lengths))]
        else:
            if verbose:
                print(f"\r gathering f-poly factors for projection: {projection_number}.", end="                                   ")
            indepSet = self.indepSets[projection_number]
            # 0 dim slice
            X = self.ring.variables
            U = tuple(var for is_indep, var in zip(indepSet, X) if is_indep)
            smallest_fpoly_factors = []
            try:
                smallest_fpoly_factors += self.extension_contraction_fpoly(U, "lp")
            except TimeoutError:
                smallest_fpoly_factors += ['- TIMEDOUT - probably very very long ' * 80]
            except Exception as e:
                if verbose:
                    print(f"An exception occurred: {e}")
                smallest_fpoly_factors += ['- TIMEDOUT - probably very very long ' * 80]
        if smallest_fpoly_factors == ['- TIMEDOUT - probably very very long ' * 80]:
            raise Inconclusive("Timedout on fpoly factors gathering.")
        smallest_fpoly_factors = sorted(smallest_fpoly_factors, key=lambda x: len(x))
        if verbose:
            print(f"\r smallest f poly factors ({len(smallest_fpoly_factors)}): {smallest_fpoly_factors}", end="                    \n")
        # check that the dimensionality drops when adding each of these factors separately (and hence drops for <ideal, f^s>)
        syngular.TIMEOUT = timeout_dim
        self.codim  # just cache it
        for i, factor in enumerate(smallest_fpoly_factors):
            if verbose:
                print(f"\r at factor {i}: {factor}.", end="                                       \n")
            X = deepcopy(self)
            X.generators += [factor]
            X.delete_cached_properties()

            # Experimental - Assumes codim w/ deg bound <= true codim.
            # Helps termiante the prime test early, IF the result is True.
            for deg_bound in [4, 6, 8, 10, 12, 14, 16, 18, 20, 22, 24] * iterated_degbound_computation:
                syngular.DEGBOUND = deg_bound
                X.delete_cached_properties()
                if self.codim < X.codim:
                    if verbose:
                        print(f"deg_bound {deg_bound} computation was conclusive.", end="\n")
                    break
                else:
                    if verbose:
                        print(f"deg_bound {deg_bound} computation was inconclusive.", end="\n")

            else:  # loop completed without encountering a break
                if verbose:
                    print("deg_bound reset to zero. Performing full computation.", end="\n")
                syngular.DEGBOUND = 0
                # if self.indepSet.count(0) >= X.indepSet.count(0):   # deprecated, equivalent to next line.
                if self.codim >= X.codim:
                    return False
        syngular.DEGBOUND = 0  # Reset it to zero
        return True

    def test_primality(self, *args, **kwargs):
        """Tests if an ideal is prime. Not to be confused with primary."""
        # eventually check assumptions and decide which test to use
        return self.primeTestDLP(*args, **kwargs)


class Inconclusive(Exception):
    pass
