import sympy

from syngular import Ring


def test_instantiation():
    ring0 = Ring('0', sympy.symbols(('x1', 'x2')), 'dp')
    ring1 = Ring('0', ('x1', 'x2'), 'dp')
    ring2 = Ring('0', ('x', 'y'), 'dp')
    assert ring0 == ring1
    assert not ring0 == ring2
