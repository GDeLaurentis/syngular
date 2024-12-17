TIMEOUT = 60  # seconds  # noqa
DEGBOUND = 0  # 0 = no-bound  # noqa
DEGBOUNDs = [4, 6, 8, 10, 12, 14, 16, 18, 20, 22, 24]  # noqa

from .ideal import Ideal  # noqa
from .ring import Ring  # noqa
from .qring import QuotientRing, QRing  # noqa
from .tools import SingularException  # noqa
from .field import Field  # noqa
from .polynomial import Polynomial, Monomial  # noqa
