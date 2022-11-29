import pytest
import sympy
import numpy

from syngular import Ideal, Ring, SingularException
from syngular.ideal import monomial_to_exponents, reduce


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


def test_radical():
    I = Ideal(Ring('0', ('x1', 'x2'), 'dp'), ['x1^2*x2'])
    J = Ideal(Ring('0', ('x1', 'x2'), 'dp'), ['x2'])
    L = Ideal(Ring('0', ('x1', 'x2'), 'dp'), ['x1'])
    assert I.radical == J & L


def test_primary_intersection():
    I = Ideal(Ring('0', ('x1', 'x2'), 'dp'), ['x1^2*x2'])
    J = Ideal(Ring('0', ('x1', 'x2'), 'dp'), ['x2'])
    K = Ideal(Ring('0', ('x1', 'x2'), 'dp'), ['x1^2'])
    assert I == J & K
    assert I == Ideal.intersection(*(J, I, K))
    assert I == Ideal.intersection(*(I, ))


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


def test_indepSet():
    I = Ideal(Ring('0', ('x1', 'x2'), 'dp'), ['x1^2*x2'])
    assert I.indepSet in [(1, 0), (0, 1)]


def test_indepSets():
    I = Ideal(Ring('0', ('x1', 'x2'), 'dp'), ['x1^2*x2'])
    assert I.indepSets == [(1, 0), (0, 1)]


def test_hash_and_set():
    I = Ideal(Ring('0', ('x1', 'x2'), 'dp'), ['x1^2*x2'])
    J = Ideal(Ring('0', ('x1', 'x2'), 'dp'), ['x2'])
    assert {I, I} == {I}
    assert {I, I, J} == {I, J, J} == {I, J}


def test_ideal_to_qring_and_back():
    ring = Ring('0', ('x1', 'x2'), 'dp')
    I = Ideal(ring, ['x1', 'x2'])
    J = Ideal(ring, ['x1'])
    I.to_qring(J)
    assert I.groebner_basis == I.minbase == ['x2']
    I.to_full_ring()
    assert set(I.groebner_basis) == set(I.minbase) == set(['x1', 'x2'])


def test_lead_gb_monomials():
    R = Ring('0', ('x', 'y', ), 'dp')
    I = Ideal(R, ['x^2', 'y'])
    assert I.leadGBmonomials == ['y', 'x2']


def test_ideal_reduce_poly():
    ring = Ring('0', ('x1', 'x2'), 'dp')
    I = Ideal(ring, ['x1'])
    poly = 'x1+x2'
    assert I.reduce(poly) == 'x2' == reduce(poly, I)


def test_ideal_reduce_ideal():
    ring = Ring('0', ('x1', 'x2'), 'dp')
    I = Ideal(ring, ['x1', 'x2'])
    J = Ideal(ring, ['x1'])
    K = J.reduce(I)
    assert K == Ideal(ring, ['x2'])
    L = I.reduce(J)
    assert L == Ideal(ring, ['0'])


def test_ideal_dim():
    I = Ideal(Ring('0', ('x1', 'x2'), 'dp'), [])
    assert I.dim == 2
    I = Ideal(Ring('0', ('x1', 'x2'), 'dp'), ['x1*x2'])
    assert I.dim == 1


def test_ideal_dims_and_codims():
    I = Ideal(Ring('0', ('x1', 'x2', 'x3'), 'dp'), ['(x3+1)*x1', '(x2+1)*x1'])
    assert I.dim == 2
    assert I.dims == {1, 2}
    assert I.codim == 1
    assert I.codims == {1, 2}


def test_monomial_to_exponents():
    r = Ring('0', ('x1', 'x2', 'x3'), 'dp')
    monomial = 'x1*x2^4*x3^123'
    assert numpy.all(monomial_to_exponents(r.variables, monomial) == numpy.array([1, 4, 123]))


def test_monomial_to_exponents_Singular_notation():
    r = Ring('0', ('x', 'y', 'z'), 'dp')
    monomial = 'x*y4*z123'
    assert numpy.all(monomial_to_exponents(r.variables, monomial) == numpy.array([1, 4, 123]))
