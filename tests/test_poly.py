from syngular import Monomial, Polynomial, Field

from fractions import Fraction

from pyadic import ModP


coeffs_and_monomials = [
    (Fraction(-215, 1), Monomial("wb^1*z^9*zb^7*w^5*X^10")),
    (Fraction(1156, 1), Monomial("wb^3*z^3*zb^2*w^11*X^4")),
    (Fraction(-125, 2), Monomial("z^4*zb^6*w^8*X^6")),
    (Fraction(-9857, 2), Monomial("wb^1*z^6*zb^6*w^5*X^8")),
    (Fraction(-954, 1), Monomial("wb^3*z^9*zb^6*w^4*X^9")),
    (Fraction(2, 1), Monomial("z^9*zb^8*X^9")),
    (Fraction(-4485, 2), Monomial("wb^4*z^2*zb^5*w^9*X^5")),
    (Fraction(216, 1), Monomial("wb^4*z^8*zb^6*w^2*X^11")),
    (Fraction(-7, 2), Monomial("wb^5*z^9*zb^1*w^9*X^6")),
    (Fraction(625, 2), Monomial("wb^4*z^2*zb^7*w^8*X^6"))
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
