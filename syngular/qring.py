import syngular
import mutableint

from .ring import Ring


class QuotientRing(Ring):

    def __init__(self, ring, ideal):
        super().__init__(ring.field, ring.variables, ring.ordering)
        self.ideal = ideal

    def __hash__(self):
        return hash(str(self))

    def __str__(self):
        string = f"{super().__str__()};\nideal i_ = {self.ideal};\nqring q = std(i_)"
        if syngular.DEGBOUND != mutableint.MutableInt(0):
            string += f";\ndegBound = {syngular.DEGBOUND};\noption()"
        return string

    def __repr__(self):
        return str(self)


QRing = QuotientRing
