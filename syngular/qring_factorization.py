import warnings
import time
import re
import sympy

from .field import Field, Q, Fp


class QuotientRingFactorization:

    def factorize(self, polynomial, as_list=False, timeout="default", padic_valuation=False, verbose=False):
        """Factorize in the quotient ring using a suitable polynomial localization, or in the ambient ring if none is provided."""
        from .ideal import Ideal
        from .point import RingPoint

        localization = getattr(self, "localization", None)

        if localization is None:
            warnings.warn("No quotient-ring factorization method is configured; falling back to the ambient polynomial ring.", stacklevel=2)
            factorization = self.poly_ring.factorize(polynomial, timeout=timeout)
            if as_list:
                return factorization
            return sympy.prod(sympy.sympify(factor.replace("^", "**")) ** multiplicity for factor, multiplicity in factorization)

        timings = {}

        if verbose > 1:
            print("localizing...")
        start = time.perf_counter()
        numerator, denominator_power = localization(polynomial, timeout=timeout)
        field = Q if self.field in (0, "0") else Fp(int(self.field))
        unit = field(1)
        element_factors = localization.element_factors
        element_powers = {factor: -denominator_power * power for factor, power in element_factors}
        timings["localization"] = time.perf_counter() - start

        if verbose > 1:
            print("factorizing...")
        start = time.perf_counter()
        factorization = self.poly_ring.factorize(numerator, timeout=timeout)
        timings["factorize"] = time.perf_counter() - start

        start = time.perf_counter()
        factors = []

        ring_symbols = set(map(str, self.poly_ring.variables))

        if verbose > 1:
            print("saturating...")
        if padic_valuation:
            Qp = Field("padic", 2 ** 31 - 1, 5)
            valuation_points = {}
            for i, (element_factor, _) in enumerate(element_factors):
                ideal = Ideal(self, (element_factor,))
                point = RingPoint(self, Qp, val=ideal.point_on_variety(field=Qp, directions=(element_factor,), valuations=(1,), seed=i))
                assert point(element_factor).n == 1
                valuation_points[element_factor] = point

        for factor, multiplicity in factorization:
            symbols = set(re.findall(r"[A-Za-z_][A-Za-z0-9_]*", factor))

            if not symbols & ring_symbols:
                factor_num, factor_den = Q(factor).as_integer_ratio()
                unit *= Q(factor_num, factor_den) ** multiplicity
                continue

            for element_factor, _ in element_factors:
                if padic_valuation:
                    valuation = valuation_points[element_factor](factor).n
                else:
                    ideal = Ideal(self, (factor,))
                    valuation = ideal.saturation_index(element_factor)
                element_powers[element_factor] += multiplicity * valuation

                if valuation != 0:
                    factor = str(self.exact_divide(factor, element_factor ** valuation))

            symbols = set(re.findall(r"[A-Za-z_][A-Za-z0-9_]*", factor))
            if symbols & ring_symbols:
                factors.append((factor, multiplicity))
            else:
                factor_num, factor_den = Q(factor).as_integer_ratio()
                unit *= Q(factor_num, factor_den) ** multiplicity

        timings["saturation"] = time.perf_counter() - start

        if any(power < 0 for power in element_powers.values()):
            raise RuntimeError(f"Negative localization valuations: {element_powers}")

        if verbose > 1:
            print("checking...")
        start = time.perf_counter()
        point = self.random_point(Field("finite field", 2 ** 31 - 1 if field is Q else field.characteristic, 1), seed=0)

        evaluated_factorization = point.field(str(unit))
        for element_factor, power in element_powers.items():
            evaluated_factorization *= point(str(element_factor)) ** power
        for factor, multiplicity in factors:
            evaluated_factorization *= point(factor) ** multiplicity

        if point(polynomial) != evaluated_factorization:
            raise RuntimeError(
                f"Factorization verification failed.\n"
                f"unit: {unit}\n"
                f"element_powers: {element_powers}\n"
                f"factors: {factors}\n"
                f"is not a valid factorization for\n{polynomial}"
            )

        timings["check"] = time.perf_counter() - start

        if verbose:
            print("polynomial:", polynomial)
            print("numerator:", numerator)
            print("factorization:", factorization)
            print("denominator_power:", denominator_power)
            print("unit:", unit)
            print("element_powers:", element_powers)
            print("factors:", factors)

            for name, timing in timings.items():
                print(f"{name:17s}: {timing:.3f} s")

        if as_list:
            return [(str(unit), 1), ] + [(str(element_factor), power) for element_factor, power in element_powers.items() if power > 0] + \
                [(str(factor), multiplicity) for factor, multiplicity in factors]

        assembled_factorization = sympy.sympify(str(unit))
        for element_factor, power in element_powers.items():
            assembled_factorization *= element_factor ** power
        for factor, multiplicity in factors:
            assembled_factorization *= sympy.sympify(factor.replace("^", "**")) ** multiplicity

        return assembled_factorization
