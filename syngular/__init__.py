import mutableint

TIMEOUT = mutableint.MutableInt(60)  # seconds  # noqa
DEGBOUND = mutableint.MutableInt(0)  # 0 = no-bound  # noqa

from .ideal import Ideal  # noqa
from .ring import Ring  # noqa
from .qring import QuotientRing  # noqa
from .tools import SingularException  # noqa
