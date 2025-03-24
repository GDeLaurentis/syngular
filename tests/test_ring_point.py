import sympy

from pyadic import ModP, PAdic
from syngular import Field, Ring, RingPoint


Qp = Field("padic", 2 ** 31 - 1, 15)
Fp = Field("finite field", 2 ** 31 - 1, 1)


def test_Fp_ring_point():
    ring = Ring('0', ('z', 'zb', 'w', 'wb', 'X'), 'dp')
    oPoint = RingPoint(ring, Fp)
    assert isinstance(oPoint('z*zb+X'), ModP)


def test_Qp_ring_point():
    ring = Ring('0', ('z', 'zb', 'w', 'wb', 'X'), 'dp')
    oPoint = RingPoint(ring, Qp)
    assert isinstance(oPoint('z*zb+X'), PAdic)


def test_Fp_ring_slice():
    ring = Ring('0', ('z', 'zb', 'w', 'wb', 'X'), 'dp')
    oSlice = RingPoint(ring, Fp, seed=0)
    oSlice.univariate_slice()
    assert isinstance(oSlice("z"), sympy.Expr)
    assert oSlice("z").as_poly().degree() == 1
