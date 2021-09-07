from point import Point2D, Point3D


class ModelPoint2D(Point2D):
    def __init__(self, x, y):
        super(ModelPoint2D, self).__init__(x, y)
        self._s = None
        self._t = None
        self._h = None
        self._v = None
        self._cell_x = None
        self._cell_y = None

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

    @property
    def cell_x(self):
        return self._cell_x

    @cell_x.setter
    def cell_x(self, index):
        self._cell_x = index

    @property
    def cell_y(self):
        return self._cell_y

    @cell_y.setter
    def cell_y(self, index):
        self._cell_y = index


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
