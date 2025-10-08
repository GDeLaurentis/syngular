import numpy


class RingPoints(list):

    def __init__(self, *args):
        list.__init__(self, *args)

    def __getattr__(self, name):
        def method(*args, **kwargs):
            return RingPoints([getattr(point, name)(*args, **kwargs) for point in self])
        return method

    def __call__(self, string_expr):
        return numpy.array([point(string_expr) for point in self])

    def __getitem__(self, some_slice):
        if isinstance(some_slice, int):
            return super().__getitem__(some_slice)
        else:
            return RingPoints(super().__getitem__(some_slice))

    @property
    def field(self):
        return self[0].field

    def __hash__(self):
        return hash(tuple(self))
