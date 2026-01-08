import sympy

from pyadic import ModP, PAdic
from syngular import Field, Ring, RingPoint


Qp = Field("padic", 2 ** 31 - 1, 15)
Fp = Field("finite field", 2 ** 31 - 1, 1)


def test_Fp_ring_point():
    ring = Ring('0', ('z', 'zb', 'w', 'wb', 'X'), 'dp')
    oPoint = RingPoint(ring, Fp)
    assert isinstance(oPoint('z*zb+X'), ModP)


def test_ring_point_parser():
    ring = Ring('0', ('z', 'zb', 'w', 'wb', 'X'), 'dp')
    oPoint = RingPoint(ring, Fp)
    assert oPoint('(w²+w·X·zb+X·z·zb(X·zb+1))') == oPoint('(w**2+w*X*zb+X*z*zb*(X*zb+1))')


def test_ring_point_parser_white_spaces():
    ring = Ring('0', ('z', 'zb', 'w', 'wb', 'X'), 'dp')
    oPoint = RingPoint(ring, Fp)
    assert oPoint("-( 6 ( -1 + w) ( +1 + z) ( wb - zb + wb zb - X z zb))") == oPoint("-(6*(-1+w)*(+1+z)*(wb-zb+wb*zb-X*z*zb))")


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
