import pytest
import sympy

from syngular import Ideal, Ring, QuotientRing, SingularException


def test_ideal_instantiation():
    I = Ideal(Ring('0', ('x1', 'x2'), 'dp'), ['x1*x2'])
    assert I is not None


def test_invalid_ideal_instantiation():
    with pytest.raises(SingularException):
        Ideal(Ring('0', ('x1', 'x2'), 'dp'), ['x*y'])


def test_ideal_quotient():
    I = Ideal(Ring('0', ('x1', 'x2'), 'dp'), ['x1*x2'])
    J = Ideal(Ring('0', ('x1', 'x2'), 'dp'), ['x1^2*x2'])
    K = Ideal(Ring('0', ('x1', 'x2'), 'dp'), ['x1'])
    R = Ideal(Ring('0', ('x1', 'x2'), 'dp'), ['1'])
    assert R / I == R
    assert I / R == I
    assert J / I == K


def test_ideal_mul_and_power():
    I = Ideal(Ring('0', ('x1', 'x2'), 'dp'), ['x1', 'x2'])
    J = Ideal(Ring('0', ('x1', 'x2'), 'dp'), ['x1^2', 'x1*x2', 'x2^2'])
    R = Ideal(Ring('0', ('x1', 'x2'), 'dp'), ['1'])
    assert I ** 2 == I * I == J
    assert I ** 0 == J ** 0 == R


def test_primary_decomposition():
    I = Ideal(Ring('0', ('x1', 'x2'), 'dp'), ['x1^2*x2'])
    J = Ideal(Ring('0', ('x1', 'x2'), 'dp'), ['x2'])
    K = Ideal(Ring('0', ('x1', 'x2'), 'dp'), ['x1^2'])
    L = Ideal(Ring('0', ('x1', 'x2'), 'dp'), ['x1'])
    assert I.primary_decomposition == [(J, J), (K, L)]


def test_primary_intersection():
    I = Ideal(Ring('0', ('x1', 'x2'), 'dp'), ['x1^2*x2'])
    J = Ideal(Ring('0', ('x1', 'x2'), 'dp'), ['x2'])
    K = Ideal(Ring('0', ('x1', 'x2'), 'dp'), ['x1^2'])
    assert I == J & K


def test_eliminate():
    I = Ideal(Ring('0', ('x1', 'x2'), 'dp'), ['x1', 'x2'])
    J = Ideal(Ring('0', ('x1', ), 'dp'), ['x1'])
    K = Ideal(Ring('0', ('x2', ), 'dp'), ['x2'])
    assert I.eliminate(range(1, 1)) == K
    assert I.eliminate(range(2, 2)) == J


def test_add():
    I = Ideal(Ring('0', ('x1', 'x2'), 'dp'), ['x1', 'x2'])
    J = Ideal(Ring('0', ('x1', 'x2'), 'dp'), ['x1'])
    K = Ideal(Ring('0', ('x1', 'x2'), 'dp'), ['x2'])
    assert I == J + K


def test_poly_contains():
    K = Ideal(Ring('0', ('x1', 'x2'), 'dp'), ['x1^2'])
    L = Ideal(Ring('0', ('x1', 'x2'), 'dp'), ['x1'])
    assert 'x1' not in K
    assert 'x1' in L
    assert sympy.symbols('x1') in L


def test_ideal_contains():
    K = Ideal(Ring('0', ('x1', 'x2'), 'dp'), ['x1^2'])
    L = Ideal(Ring('0', ('x1', 'x2'), 'dp'), ['x1'])
    assert K in L
    assert L not in K


def test_indepSets():
    I = Ideal(Ring('0', ('x1', 'x2'), 'dp'), ['x1^2*x2'])
    assert I.indepSets == [(1, 0), (0, 1)]


def test_hash_and_set():
    I = Ideal(Ring('0', ('x1', 'x2'), 'dp'), ['x1^2*x2'])
    J = Ideal(Ring('0', ('x1', 'x2'), 'dp'), ['x2'])
    assert {I, I} == {I}
    assert {I, I, J} == {I, J, J} == {I, J}


def test_ideal_over_qring():
    ring = Ring('0', ('x1', 'x2'), 'dp')
    I = Ideal(ring, ['x1', 'x2'])
    J = Ideal(ring, ['x1'])
    qring = QuotientRing(ring, J)
    I.ring = qring
    assert I.groebner_basis == I.minbase == ['x2']


def test_ideal_reduce():
    ring = Ring('0', ('x1', 'x2'), 'dp')
    I = Ideal(ring, ['x1'])
    poly = 'x1+x2'
    remainder = I.reduce(poly)
    assert remainder == 'x2'
