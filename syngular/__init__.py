TIMEOUT = 60  # seconds  # noqa
DEGBOUND = 0  # noqa, 0 = no-bound  
DEGBOUNDs = [4, 6, 8, 10, 12, 14, 16, 18, 20, 22, 24]  # noqa
DEBUG = False  # noqa
FORCECDOTS = False  # noqa
CDOTCHAR = '·'  # noqa, '·' or '*' or ' '
UNICODEPOWERS = True  # noqa, True or False
NORMALIZE_POWERS_PATTERNS = ()  # eg. re.compile(r"(mt)(\d+)") will force mt2 to be treated as mt^2
USE_ELLIPSIS_FOR_PRINT = False  # noqa, toggles ellipsis in for str. Use locally for prints only.

from .version import __version__  # noqa
from .ideal import Ideal  # noqa
from .ring import Ring  # noqa
from .qring import QuotientRing, QRing  # noqa
from .tools import SingularException, Singular_version, flatten  # noqa
from .field import Field, Q, Qi  # noqa
from .polynomial import Polynomial, Monomial  # noqa
from .settings import TemporarySetting  # noqa
from .point import RingPoint  # noqa
from .points import RingPoints  # noqa
