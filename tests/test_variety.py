import numpy
import sympy

from copy import copy

from syngular import Field, Ring, QRing, Ideal, Polynomial


def test_finite_field_variety_point():
    Fp = Field("finite field", 2 ** 31 - 1, 1)
    I = Ideal(Ring('0', ('z', 'zb', 'w', 'wb', 'X'), 'dp'), ['zb^2*wb*X^2-zb*w*wb*X-zb^2*X^2-zb*wb*X^2+zb*wb*X+zb*X^2-w*wb+wb',
                                                             'z*w*wb*X+zb*w*wb*X+z*w^2+z*w*wb-z*w*X-zb*w*X-z*wb*X-w*wb*X-2*z*w-w^2-z*wb+z*X+w*X+z+w',
                                                             'z*zb*X-z*w+z+w'])
    point_dict = I.point_on_variety(Fp)
    assert numpy.all(numpy.array(I.generators_eval(**point_dict)) == 0)


def test_padic_variety_point():
    Qp = Field("padic", 2 ** 31 - 1, 10)
    I = Ideal(Ring('0', ('z', 'zb', 'w', 'wb', 'X'), 'dp'), ['zb^2*wb*X^2-zb*w*wb*X-zb^2*X^2-zb*wb*X^2+zb*wb*X+zb*X^2-w*wb+wb',
                                                             'z*w*wb*X+zb*w*wb*X+z*w^2+z*w*wb-z*w*X-zb*w*X-z*wb*X-w*wb*X-2*z*w-w^2-z*wb+z*X+w*X+z+w',
                                                             'z*zb*X-z*w+z+w'])
    point_dict = I.point_on_variety(Qp, directions=(
        'X**2*wb*zb**2/123 - X**2*wb*zb/123 - X**2*zb**2/123 + X**2*zb/123 - X*w*wb*zb/123 + X*wb*zb/123 - w*wb/123 + wb/123',
        'X*z*zb/456 - w*z/456 + w/456 + z/456'),
        valuations=(1, 2)
    )
    assert numpy.all(numpy.array([entry.n for entry in I.generators_eval(**point_dict)]) >= 1)
    assert eval('(zb^2*wb*X^2-zb*w*wb*X-zb^2*X^2-zb*wb*X^2+zb*wb*X+zb*X^2-w*wb+wb)/123'.replace("^", "**"), point_dict).n == 1
    assert eval('(z*zb*X-z*w+z+w)/456'.replace("^", "**"), point_dict).n == 2


def test_complex_variety_point():
    C = Field("mpc", 0, 300)
    I = Ideal(Ring('0', ('z', 'zb', 'w', 'wb', 'X'), 'dp'), ['zb^2*wb*X^2-zb*w*wb*X-zb^2*X^2-zb*wb*X^2+zb*wb*X+zb*X^2-w*wb+wb',
                                                             'z*w*wb*X+zb*w*wb*X+z*w^2+z*w*wb-z*w*X-zb*w*X-z*wb*X-w*wb*X-2*z*w-w^2-z*wb+z*X+w*X+z+w',
                                                             'z*zb*X-z*w+z+w'])
    point_dict = I.point_on_variety(C, valuations=(10 ** -30, 10 ** -30))
    assert numpy.all(numpy.isclose(numpy.array(I.generators_eval(**point_dict)).astype(complex), 0))


def test_complex_variety_point_with_base_point_and_unequal_valuations():
    C = Field("mpc", 0, 300)
    I = Ideal(Ring('0', ('z', 'zb', 'w', 'wb', 'X'), 'dp'), ['zb^2*wb*X^2-zb*w*wb*X-zb^2*X^2-zb*wb*X^2+zb*wb*X+zb*X^2-w*wb+wb',
                                                             'z*w*wb*X+zb*w*wb*X+z*w^2+z*w*wb-z*w*X-zb*w*X-z*wb*X-w*wb*X-2*z*w-w^2-z*wb+z*X+w*X+z+w',
                                                             'z*zb*X-z*w+z+w'])
    base_point = {'X': 1, 'z': 2, 'zb': 3, 'w': 5, 'wb': 7}
    point_dict = I.point_on_variety(C, base_point=base_point,
                                    directions=('zb^2*wb*X^2-zb*w*wb*X-zb^2*X^2-zb*wb*X^2+zb*wb*X+zb*X^2-w*wb+wb', 'z*zb*X-z*w+z+w'),
                                    valuations=(10 ** -30, 10 ** -60))
    assert numpy.all(numpy.isclose(numpy.array(I.generators_eval(**point_dict)).astype(complex), 0))
    assert numpy.isclose(abs(complex(eval('zb^2*wb*X^2-zb*w*wb*X-zb^2*X^2-zb*wb*X^2+zb*wb*X+zb*X^2-w*wb+wb'.replace("^", "**"), copy(point_dict)))), 10 ** -30)
    assert numpy.isclose(abs(complex(eval('z*zb*X-z*w+z+w'.replace("^", "**"), copy(point_dict)))), 10 ** -60)
    assert [point_dict[key] == base_point[key] for key in base_point.keys()].count(True) == 3


def test_padic_variety_point_in_qring():
    """This is ("Δ_14|23|56", "⟨1|2+3|4]" ) from lips at six point."""
    ring = Ring(0, ('a1', 'b1', 'c1', 'd1', 'a2', 'b2', 'c2', 'd2', 'a3', 'b3', 'c3', 'd3', 'a4', 'b4', 'c4', 'd4', 'a5', 'b5', 'c5', 'd5', 'a6', 'b6', 'c6', 'd6'), 'dp')
    ideal = Ideal(ring, ['b1*d1 + b2*d2 + b3*d3 + b4*d4 + b5*d5 + b6*d6', '-a1*d1 - a2*d2 - a3*d3 - a4*d4 - a5*d5 - a6*d6',
                         '-b1*c1 - b2*c2 - b3*c3 - b4*c4 - b5*c5 - b6*c6', 'a1*c1 + a2*c2 + a3*c3 + a4*c4 + a5*c5 + a6*c6'])
    qring = QRing(ring, ideal)
    directions = [
        '-((a1*c1 + a4*c4)*(b1*d1 + b4*d4) - (-a1*d1 - a4*d4)*(-b1*c1 - b4*c4))*((a2*c2 + a3*c3)*(b2*d2 + b3*d3) - (a2*d2 + a3*d3)*(b2*c2 + b3*c3)) + ((a1*c1 + a4*c4)*(b2*d2 + b3*d3)/2 + (-a1*d1 - a4*d4)*(b2*c2 + b3*c3)/2 + (a2*c2 + a3*c3)*(b1*d1 + b4*d4)/2 + (a2*d2 + a3*d3)*(-b1*c1 - b4*c4)/2)**2',  # noqa
        '-c4*(-a1*(b2*d2 + b3*d3) + b1*(a2*d2 + a3*d3)) + d4*(-a1*(b2*c2 + b3*c3) + b1*(a2*c2 + a3*c3))'
    ]
    directions = [sympy.sympify(direction) for direction in directions]
    J = Ideal(qring, directions)
    Qp = Field("padic", 2 ** 31 - 1, 11)
    point = J.point_on_variety(Qp, directions=directions, valuations=(2, 2))
    assert all([Polynomial(direction.expand(), Field("rational", 0, 0))(point, Qp).coeffs[0].n == 2 for direction in directions])
    # This is (s_123-s_234), which lives in the radical of J but not in J.
    member_of_radical = '(a1*c1 + a2*c2 + a3*c3)*(b1*d1 + b2*d2 + b3*d3) - (-a1*d1 - a2*d2 - a3*d3)*(-b1*c1 - b2*c2 - b3*c3) - (a2*c2 + a3*c3 + a4*c4)*(b2*d2 + b3*d3 + b4*d4) + (-a2*d2 - a3*d3 - a4*d4)*(-b2*c2 - b3*c3 - b4*c4)'  # noqa
    member_of_radical = sympy.sympify(member_of_radical).expand()
    assert Polynomial(member_of_radical.expand(), Field("rational", 0, 0))(point, Qp).coeffs[0].n == 1
