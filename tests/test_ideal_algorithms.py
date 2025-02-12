import pytest
import syngular

from syngular import Ideal, Ring
from syngular.ideal_algorithms import Inconclusive


def test_extension_contraction():
    I = Ideal(Ring('0', ('x1', 'x2'), 'dp'), ['x1*x2'])
    J = Ideal(Ring('0', ('x1', 'x2'), 'dp'), ['x2'])
    indepSet = I.indepSets[0]
    X = I.ring.variables
    U = tuple(var for is_indep, var in zip(indepSet, X) if is_indep)
    assert I.extension_contraction(U, ordering="lp") == (1, J)
    assert I.extension_contraction(U, ordering="dp") == (1, J)
    with pytest.raises(ValueError):
        I.extension_contraction(U, ordering="Dp") == (1, J)


def test_test_primality():
    I = Ideal(Ring('0', ('x1', 'x2'), 'dp'), ['x1*x2'])
    J = Ideal(Ring('0', ('x1', 'x2'), 'dp'), ['x2'])
    assert I.test_primality(verbose=True) is False
    assert I.test_primality(verbose=True, iterated_degbound_computation=True) is False
    assert J.test_primality(verbose=True) is True
    assert J.test_primality(verbose=True, iterated_degbound_computation=True) is True


def test_test_primality_primary_but_not_prime():
    I = Ideal(Ring('0', ('x1', 'x2'), 'dp'), ['x1^2'])
    assert I.test_primality(verbose=True) is False
    assert I.test_primality(verbose=True, iterated_degbound_computation=True) is False
    assert I.test_primality(verbose=True, astuple=True) == (True, False)
    assert I.test_primality(verbose=True, iterated_degbound_computation=True, astuple=True) == (True, False)


def test_test_primality_with_degree_bound():
    # from: LipsIdeal(6, ('⟨1|2+3|1]', '⟨2|1+4|2]'))
    r = Ring(0, ('d6', 'c6', 'b6', 'a6', 'd5', 'c5', 'b5', 'a5', 'd4', 'c4', 'b4', 'a4', 'd3', 'c3', 'b3', 'a3', 'd2', 'c2', 'b2', 'a2', 'd1', 'c1', 'b1', 'a1'), 'dp')
    generators = ['4*a1*b2*c1*d2 - 4*a1*b2*c2*d1 + 4*a1*b3*c1*d3 - 4*a1*b3*c3*d1 - 4*a2*b1*c1*d2 + 4*a2*b1*c2*d1 - 4*a3*b1*c1*d3 + 4*a3*b1*c3*d1',
                  '4*a1*b2*c1*d2 - 4*a1*b2*c2*d1 - 4*a2*b1*c1*d2 + 4*a2*b1*c2*d1 + 4*a2*b4*c2*d4 - 4*a2*b4*c4*d2 - 4*a4*b2*c2*d4 + 4*a4*b2*c4*d2',
                  'b1*d1 + b2*d2 + b3*d3 + b4*d4 + b5*d5 + b6*d6',
                  '-a1*d1 - a2*d2 - a3*d3 - a4*d4 - a5*d5 - a6*d6',
                  '-b1*c1 - b2*c2 - b3*c3 - b4*c4 - b5*c5 - b6*c6',
                  'a1*c1 + a2*c2 + a3*c3 + a4*c4 + a5*c5 + a6*c6']
    I = Ideal(r, generators)
    syngular.DEGBOUNDs = [4, 6, ]
    assert I.test_primality(verbose=True, iterated_degbound_computation=True, projection_number=192)
