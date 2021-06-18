from .ring import Ring


class QuotientRing(Ring):

    def __init__(self, ring, ideal):
        super().__init__(ring.field, ring.variables, ring.ordering)
        self.ideal = ideal

    def __hash__(self):
        return hash(str(self))

    def __str__(self):
        return f"{super().__str__()};\nideal i_ = {self.ideal};\nqring q = std(i_)"

    def __repr__(self):
        return str(self)
