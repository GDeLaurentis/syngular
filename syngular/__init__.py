TIMEOUT = 60  # seconds  # noqa
DEGBOUND = 0  # 0 = no-bound  # noqa
DEGBOUNDs = [4, 6, 8, 10, 12, 14, 16, 18, 20, 22, 24]  # noqa
DEBUG = False  # noqa

from .version import __version__  # noqa
from .ideal import Ideal  # noqa
from .ring import Ring  # noqa
from .qring import QuotientRing, QRing  # noqa
from .tools import SingularException, Singular_version, flatten  # noqa
from .field import Field  # noqa
from .polynomial import Polynomial, Monomial  # noqa
from .settings import TemporarySetting  # noqa
from .point import RingPoint  # noqa
