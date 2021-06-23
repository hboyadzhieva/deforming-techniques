import math

from point import Point2D, Point3D


class ModelPoint2D(Point2D):
    def __init__(self, x, y):
        super(ModelPoint2D, self).__init__(x, y)
        self._s = None
        self._t = None
        self._h = None
        self._v = None

    @property
    def s(self):
        return self._s

    @s.setter
    def s(self, s):
        self._s = s

    @property
    def t(self):
        return self._t

    @t.setter
    def t(self, t):
        self._t = t

    @property
    def h(self):
        return self._h

    @h.setter
    def h(self, h):
        self._h = h

    @property
    def v(self):
        return self._v

    @v.setter
    def v(self, v):
        self._v = v

    def right(self):
        if self.s is not None:
            return math.ceil(self.s)

    def left(self):
        if self.s is not None:
            return math.ceil(self.s) - 1

    def top(self):
        if self.t is not None:
            return math.ceil(self.t)

    def bottom(self):
        if self.t is not None:
            return math.ceil(self.t) - 1


class ModelPoint3D(Point3D):
    def __init__(self, x, y, z):
        super(ModelPoint3D, self).__init__(x, y, z)
        self._s = None
        self._t = None
        self._u = None

    @property
    def s(self):
        return self._s

    @s.setter
    def s(self, s):
        self._s = s

    @property
    def t(self):
        return self._t

    @t.setter
    def t(self, t):
        self._t = t

    @property
    def u(self):
        return self._u

    @u.setter
    def u(self, u):
        self._u = u
