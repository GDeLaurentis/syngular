from syngular import Field, Ring, RingPoint, RingPoints


Fp = Field("finite field", 2 ** 31 - 1, 1)
ring = Ring('0', ('z', 'zb', 'w', 'wb', 'X'), 'dp')


def test_ring_points_method_forwarding_and_call():
    lPs = [RingPoint(ring, field=Fp, seed=seed) for seed in range(10)]
    lPs = RingPoints(lPs)
    lPs.univariate_slice()
    assert all('t' in str(entry) for entry in lPs('X'))


def test_ring_points_indexing():
    lPs = [RingPoint(ring, field=Fp, seed=seed) for seed in range(10)]
    lPs = RingPoints(lPs)
    assert [lPs[0]('X')] == lPs[:1]('X')


def test_ring_points_field():
    lPs = [RingPoint(ring, field=Fp, seed=seed) for seed in range(10)]
    lPs = RingPoints(lPs)
    assert lPs.field == lPs[0].field
