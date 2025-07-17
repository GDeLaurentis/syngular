import numpy


class RingPoints(list):

    def __init__(self, *args):
        list.__init__(self, *args)

    def __getattr__(self, name):
        def method(*args, **kwargs):
            return RingPoints([getattr(oPs, name)(*args, **kwargs) for oPs in self])
        return method

    def __call__(self, string_expr):
        return numpy.array([oPs(string_expr) for oPs in self])

    def __getitem__(self, some_slice):
        if isinstance(some_slice, int):
            return super().__getitem__(some_slice)
        else:
            return RingPoints(super().__getitem__(some_slice))

    @property
    def field(self):
        return self[0].field
