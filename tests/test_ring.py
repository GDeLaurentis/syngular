import numpy
import sympy
import pytest

from syngular import Ideal, Ring, QRing, Field, Polynomial


def test_instantiation():
    ring0 = Ring('0', sympy.symbols(('x1', 'x2')), 'dp')
    ring1 = Ring('0', ('x1', 'x2'), 'dp')
    ring2 = Ring('0', ('x', 'y'), 'dp')
    assert ring0 == ring1
    assert not ring0 == ring2


def test_trivial_instantiation():
    Ring('0', ('x', ), 'dp')


def test_qring_and_point_in_ring():
    m, D = 6, 7
    η = numpy.diag([1, ] + [-1, ] * (D - 1))
    momenta = numpy.array([numpy.array([sympy.symbols(f"p{i}_{j}") for j in range(D)]) for i in range(1, m + 1)])
    on_shell_relations = [momentum @ η @ momentum for momentum in momenta]
    momentum_conservation = sum(momenta).tolist()
    r = Ring('0', tuple(momenta.flatten().tolist()), 'dp')
    i = Ideal(r, on_shell_relations + momentum_conservation)
    q = QRing(r, i)
    field = Field('finite field', 2 ** 31 - 1, 1)
    point = q.random_point(field=field)
    assert all([Polynomial(generator, field)(point) == 0 for generator in i.generators])


def test_exact_divide():
    ring = Ring(0, ("x", "y"), "dp")
    quotient = ring.exact_divide("(x+y)*(x-y)", "x+y")
    x, y = sympy.symbols("x y")
    assert sympy.expand(sympy.sympify(quotient.replace("^", "**")) - (x - y)) == 0


def test_exact_divide_not_divisible():
    ring = Ring(0, ("x", "y"), "dp")
    with pytest.raises(ArithmeticError):
        ring.exact_divide("x^2+y", "x+y")


def test_factorize():
    ring = Ring(0, ("x", "y"), "dp")
    factors = ring.factorize("(x-y)^2*(x+2*y)^3")
    x, y = sympy.symbols("x y")
    reconstructed = sympy.prod(sympy.sympify(factor.replace("^", "**")) ** power for factor, power in factors)
    assert sympy.expand(reconstructed - (x - y) ** 2 * (x + 2 * y) ** 3) == 0
    assert sorted(power for _, power in factors) == [1, 2, 3]


def test_singular_factorize_in_qring_warns_and_fails():
    ring = Ring(0, ("x", "y", "z"), "dp")
    qring = QRing(ring, Ideal(ring, ("x*y-z",)))
    with pytest.warns(UserWarning, match="ambient polynomial ring"):
        assert qring.factorize("z") != qring.factorize("x*y")


def test_zero_and_unit_ideal_in_quotient_ring():
    ring = Ring(0, ("x", "y", "z"), "dp")
    qring = QRing(ring, Ideal(ring, ("x*y-z",)))
    assert qring.zero_ideal() == Ideal(qring, ("x*y-z",))
    assert qring.unit_ideal() == Ideal(qring, ("1+x*y-z",))


def test_qring_factorize():
    ring = Ring(0, ("x", "y", "z"), "dp")
    qring = QRing(ring, Ideal(ring, ("z-x*y",)))
    qring.set_localization("x", dependent_variables=("y",), check_polynomial=True)

    x, y = sympy.symbols("x y")

    assert sympy.expand(qring.factorize("z") - x * y) == 0
    assert sympy.expand(qring.factorize("z", padic_valuation=True) - x * y) == 0
    assert sympy.expand(qring.factorize("z*(z+x)") - x**2 * y * (y + 1)) == 0
