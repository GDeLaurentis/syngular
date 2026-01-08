from .version import __version__
from .ideal import Ideal
from .ring import Ring
from .qring import QuotientRing, QRing
from .tools import SingularException, Singular_version
from .field import Field, Q, Qi
from .polynomial import Polynomial, Monomial
from .point import RingPoint
from .points import RingPoints


TIMEOUT = 60  # seconds  # noqa
DEGBOUND = 0  # noqa, 0 = no-bound  
DEGBOUNDs = [4, 6, 8, 10, 12, 14, 16, 18, 20, 22, 24]  # noqa
DEBUG = False  # noqa
FORCECDOTS = False  # noqa
CDOTCHAR = '·'  # noqa, '·' or '*' or ' '
UNICODEPOWERS = True  # noqa, True or False
NORMALIZE_POWERS_PATTERNS = ()  # eg. re.compile(r"(mt)(\d+)") will force mt2 to be treated as mt^2
USE_ELLIPSIS_FOR_PRINT = False  # noqa, toggles ellipsis in for str. Use locally for prints only.

__all__ = [
    "__version__",
    "Ideal",
    "Ring",
    "QuotientRing",
    "QRing",
    "SingularException",
    "Singular_version",
    "Field",
    "Q",
    "Qi",
    "Polynomial",
    "Monomial",
    "TemporarySetting",
    "RingPoint",
    "RingPoints",
]


# Back-compatibility - to be removed

import warnings  # noqa


def __getattr__(name):
    if name in {"flatten", }:
        warnings.warn(
            f"syngular.{name} is deprecated and will be removed in a future release; "
            f"use pycoretools.iterables.{name} (or: from pycoretools import {name}).",
            FutureWarning,
            stacklevel=2,
        )
        from pycoretools import iterables
        return getattr(iterables, name)
    elif name in {"TemporarySetting", }:
        warnings.warn(
            f"syngular.{name} is deprecated and will be removed in a future release; "
            f"use pycoretools.context.{name} (or: from pycoretools import {name}).",
            FutureWarning,
            stacklevel=2,
        )
        from pycoretools import context
        return getattr(context, name)
    raise AttributeError(f"module {__name__!r} has no attribute {name!r}")


def __dir__():
    return sorted(list(globals().keys()) + ["flatten", "TemporarySetting"])
