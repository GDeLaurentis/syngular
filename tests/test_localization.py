import sympy
from syngular import Ring, QRing, Ideal


def test_principal_localization():
    ring = Ring(0, ("x", "y", "z"), "dp")
    qring = QRing(ring, Ideal(ring, ("x*y-z",)))
    localization = qring.principal_localization("x", dependent_variables=("y",), check_polynomial=True)

    assert localization.is_polynomial_localization

    numerator, power = localization("y+1")
    x, z = sympy.symbols("x z")

    assert power == 1
    assert sympy.expand(sympy.sympify(numerator.replace("^", "**")) - (x + z)) == 0

    numerator, power = localization("x*y")
    assert power == 0
    assert sympy.expand(sympy.sympify(numerator.replace("^", "**")) - z) == 0
