import syngular

from .tools import execute_singular_command
from .ring import Ring
from .qring_factorization import QuotientRingFactorization


class QuotientRing(QuotientRingFactorization, Ring):

    def __init__(self, ring, ideal):
        if ring != ideal.ring:
            raise ValueError("The quotient-ring ambient ring and the ideal ring must agree.")
        super().__init__(ring.field, ring.variables, ring.ordering)
        self.ideal = ideal
        self.test_valid_qring()

    @property
    def poly_ring(self):
        return self.ideal.ring

    def test_valid_qring(self):
        singular_commands = [f"ring r = {QuotientRing.__str__(self)};", "print(q);", "$"]
        execute_singular_command(singular_commands)

    def __hash__(self):
        return hash(str(self))

    def __str__(self):
        string = f"{super().__str__()};\nideal i_ = {self.ideal};\nqring q = std(i_)"
        string += ";\nshort=0"
        if syngular.DEGBOUND != 0:
            string += f";\ndegBound = {syngular.DEGBOUND}"
        return string

    def __repr__(self):
        return str(self)

    def principal_localization(self, element, dependent_variables=(), check_polynomial=False):
        """Construct a principal localization of this quotient ring."""
        from .localization import PrincipalLocalization
        return PrincipalLocalization(self, element, dependent_variables, check_polynomial)

    def set_localization(self, element, dependent_variables=(), check_polynomial=False):
        """Set an explicitly specified principal localization."""
        self._localization = self.principal_localization(element, dependent_variables, check_polynomial)
        return self._localization

    def find_localization(self, independent_variables=None, dependent_variables=None, check_polynomial=True):
        """Attempt to find and set a principal localization yielding a polynomial localization for a suitable coordinate projection."""
        from .localization import find_principal_localization
        self._localization = find_principal_localization(self, independent_variables, dependent_variables, check_polynomial)
        return self._localization

    @property
    def localization(self):
        if not hasattr(self, "_localization"):
            self._set_default_localization()
        return self._localization

    def _set_default_localization(self):
        self.find_localization()


QRing = QuotientRing
