import numpy

from copy import copy

from syngular import Field, Ring, Ideal


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
    point_dict = I.point_on_variety(Qp, directions=('(zb^2*wb*X^2-zb*w*wb*X-zb^2*X^2-zb*wb*X^2+zb*wb*X+zb*X^2-w*wb+wb)/123', '(z*zb*X-z*w+z+w)/456'), valuations=(1, 2))
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
