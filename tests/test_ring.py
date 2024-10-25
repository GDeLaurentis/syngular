import numpy
import sympy

from syngular import Ideal, Ring, QRing, Field, Polynomial


def test_instantiation():
    ring0 = Ring('0', sympy.symbols(('x1', 'x2')), 'dp')
    ring1 = Ring('0', ('x1', 'x2'), 'dp')
    ring2 = Ring('0', ('x', 'y'), 'dp')
    assert ring0 == ring1
    assert not ring0 == ring2


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