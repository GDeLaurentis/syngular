from syngular import Monomial, Polynomial, Field, Ring, RingPoint

from fractions import Fraction

from pyadic import ModP


coeffs_and_monomials = [
    (Fraction(-215, 1), Monomial(" wb^1*z^9·zb^7w^5*X^10")),
    (Fraction(1156, 1), Monomial("wb^3 z³*zb^2*w^11*X^4")),
    (Fraction(-125, 2), Monomial("z^4*zb^6*w^8*X^6")),
    (Fraction(-9857, 2), Monomial("wb^1*z^6*zb^6*w^5*X^8")),
    (Fraction(-954, 1), Monomial("wb³*z^9*zb^6·w^4*X^9")),
    (Fraction(2, 1), Monomial("z^9*zb^8*X^9")),
    (Fraction(-4485, 2), Monomial("wb⁴*z^2*zb^5*w^9X^5")),
    (Fraction(216, 1), Monomial("wb^4*z^8*zb^6*w^2*X^11")),
    (Fraction(-7, 2), Monomial("wb^5*z^9*zb^1*w^9*X^6")),
    (Fraction(625, 2), Monomial("wb^4*z^2·zb^7*w^8*X^6")),
    (Fraction(625, 2), Monomial("wb·z·zb^7w^8 X^6"))
]

point = {
    'z': ModP('1978390662 % 2147481317'),
    'zb': ModP('117802826 % 2147481317'),
    'w': ModP('962111889 % 2147481317'),
    'wb': ModP('1912527506 % 2147481317'),
    'X': ModP('747011570 % 2147481317')
}

QQ = Field("rational", 0, 0)
Fp = Field('finite field', 2147481317, 1)


def test_eq_zero():
    poly = Polynomial(coeffs_and_monomials, QQ)
    assert poly - poly == 0


def test_rstr():
    poly = Polynomial(coeffs_and_monomials, QQ)
    poly_reloaded = Polynomial(str(poly), QQ)
    assert len(poly_reloaded) == len(poly)
    assert Polynomial(poly_reloaded, QQ) == poly
    assert Polynomial(poly_reloaded, QQ) - poly == 0
    assert poly.subs(point, Fp) == poly_reloaded.subs(point, Fp)


def test_monomial_instantiation_from_tuple_of_strings():
    assert Monomial(("zb", "zb", "X²", "(wb-1)")) == Monomial("zb²X²(wb-1)")


def test_monomial_instantiation_from_tuple_of_tuples():
    assert (
        Monomial((('zb', 2), ('X', 2), ('(wb-1)', 1), ('(X z zb + 1)', 1), ('(w+z-w·z+X·z·zb)', 4))) ==
        Monomial("zb·zb·X²(wb-1)(X z zb + 1)(w+z-w·z+X·z·zb)⁴")
    )


def test_monomial_instantiation_from_pair_of_tuples():
    assert Monomial(("zb", "X", "(wb-1)"), (1, 2, 3)) == Monomial("zb·X²(wb-1)³")


def test_addition_with_field_element():
    field = Field('padic', 2 ** 31 - 19, 10)
    a, b = field.random(), field.random()
    assert (b * Polynomial([(field.one, Monomial('t')), ], field) + a ==
            a + b * Polynomial([(field.one, Monomial('t')), ], field))


def test_polynomial_with_brackets():
    string = "37/24*tr(1+2|3)⟨3|4⟩ ⟨2|3|1] + 19/24 (s_12-s_34)⟨1|2⟩·[1|2]⟨3|4⟩-23 [1|2]  ⟨2|4⟩⟨3|1+2|3|2⟩/24"
    poly = Polynomial(string, Field("rational", 0, 0))
    assert [(coeff, list(monomial.items())) for coeff, monomial in poly.coeffs_and_monomials] == [
        (Fraction(37, 24),
         [('tr(1+2|3)', 1), ('⟨3|4⟩', 1), ('⟨2|3|1]', 1)]),
        (Fraction(19, 24),
         [('(s_12-s_34)', 1), ('⟨1|2⟩', 1), ('[1|2]', 1), ('⟨3|4⟩', 1)]),
        (Fraction(-23, 24),
         [('[1|2]', 1), ('⟨2|4⟩', 1), ('⟨3|1+2|3|2⟩', 1)])]
    assert poly == Polynomial(str(poly), Field("rational", 0, 0))


def test_monomial_division():
    a = Monomial("""(zb)(X)²(wb-1)(X*z*zb+1)(w+z-w*z+X*z*zb)⁴""")
    b = Monomial("""(zb)(X)(wb-1)(zb-1)(X*z*zb+1)(w+z-w*z+X*z*zb)""")
    x = a & b
    assert a / x == Monomial("(X)(w+z-w*z+X*z*zb)³")


def test_monomial_eval_vs_ring_point_eval():
    ring = Ring('0', ('z', 'zb', 'w', 'wb', 'X'), 'dp')
    oPoint = RingPoint(ring, Field("padic", 2 ** 31 - 1, 10))
    Monomial("""zb X²(wb-1)(X z zb + 1)(w+z-w·z+X*z·zb)⁴""")(oPoint) == oPoint("zb X²(wb-1)(X z zb + 1)(w+z-w z+X·z·zb)⁴")
