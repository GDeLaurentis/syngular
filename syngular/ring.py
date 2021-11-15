import sympy
import mutableint
import syngular


class Ring(object):

    def __init__(self, field, variables, ordering):
        self.field = field
        self.variables = variables
        self.ordering = ordering

    def __hash__(self):
        return hash(str(self))

    @property
    def field(self):
        return self._field

    @field.setter
    def field(self, field):
        if isinstance(field, tuple) and isinstance(field[0], str):
            field = sympy.symbols(field)
        self._field = field

    @property
    def variables(self):
        return self._variables

    @variables.setter
    def variables(self, variables):
        if not isinstance(variables, tuple):
            raise Exception("Ring variables must be tuple")
        if len(variables) == 0:
            raise Exception("Ring must be defined over at least one variable")
        if isinstance(variables[0], str):
            variables = sympy.symbols(variables)
        self._variables = variables

    @property
    def ordering(self):
        """Monomial ordering"""
        return self._ordering

    @ordering.setter
    def ordering(self, ordering):
        if isinstance(ordering, str):
            assert ordering in ['lp', 'rp', 'dp', 'Dp']
        else:
            ordering = str(ordering).replace("'", "")
        self._ordering = ordering

    def __str__(self):
        string = ", ".join(map(str, [self.field, self.variables, self.ordering])).replace(",)", ")")
        if syngular.DEGBOUND != mutableint.MutableInt(0):
            string += f";\ndegBound = {syngular.DEGBOUND};\noption()"
        return string

    def __repr__(self):
        return str(self)

    def __eq__(self, other):
        assert isinstance(self, Ring) and isinstance(other, Ring)
        return self.field == other.field and self.variables == other.variables and self.ordering == other.ordering
