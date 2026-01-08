import functools

from pycoretools import TemporarySetting, with_cm


def with_other_cas_compatible_str(func):
    @with_cm(TemporarySetting("syngular", "FORCECDOTS", True))
    @with_cm(TemporarySetting("syngular", "CDOTCHAR", "*"))
    @with_cm(TemporarySetting("syngular", "UNICODEPOWERS", False))
    @functools.wraps(func)
    def wrapper(*args, **kwargs):
        return func(*args, **kwargs)
    return wrapper
