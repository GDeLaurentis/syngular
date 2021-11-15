import pytest

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


def test_primeTestDLP():
    I = Ideal(Ring('0', ('x1', 'x2'), 'dp'), ['x1*x2'])
    J = Ideal(Ring('0', ('x1', 'x2'), 'dp'), ['x2'])
    assert I.primeTestDLP(verbose=True) is False
    assert J.primeTestDLP(verbose=True) is True


def test_primeTestDLP_inconclusive():
    I = Ideal(Ring('0', ('x1', 'x2'), 'dp'), ['x1^2'])
    with pytest.raises(Inconclusive):
        assert I.primeTestDLP(verbose=True) is False
