import pytest
import pickle
import hashlib

from syngular import Field


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #


def test_type_casting():

    assert Field("mpc", 0, 300)(0) == 0  # what about tollerance here
    assert Field("padic", 2 ** 31 - 1, 5)(0) == 0  # and here
    assert Field("finite field", 2 ** 31 - 1, 1)(0) == 0


@pytest.mark.parametrize(
    'original', [
        Field("mpc", 0, 300),
        Field("padic", 2 ** 31 - 1, 5),
        Field("finite field", 2 ** 31 - 1, 1),
    ]
)
def test_serializable_and_hash_stable(original):
    dumped = pickle.dumps(original)
    loaded = pickle.loads(dumped)

    assert original == loaded

    hash1 = hashlib.sha256(pickle.dumps(original)).hexdigest()
    hash2 = hashlib.sha256(pickle.dumps(loaded)).hexdigest()

    assert hash1 == hash2
