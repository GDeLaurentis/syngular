import syngular

from .tools import execute_singular_command
from .ring import Ring


class QuotientRing(Ring):

    def __init__(self, ring, ideal):
        super().__init__(ring.field, ring.variables, ring.ordering)
        self.ideal = ideal
        self.test_valid_qring()

    def test_valid_qring(self):
        singular_commands = [f"ring r = {QuotientRing.__str__(self)};",
                             "print(q);"
                             "$"]
        execute_singular_command(singular_commands)

    def __hash__(self):
        return hash(str(self))

    def __str__(self):
        string = f"{super().__str__()};\nideal i_ = {self.ideal};\nqring q = std(i_)"
        if syngular.DEGBOUND != 0:
            string += f";\ndegBound = {syngular.DEGBOUND};\noption()"
        return string

    def __repr__(self):
        return str(self)


QRing = QuotientRing
