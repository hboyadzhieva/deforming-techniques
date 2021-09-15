"""
Point in n-dim space that comprises the model for deformation.
Besides the coordinates of the point, the model point also consists
of parameters, that are relative to the control grid.
"""
from point import Point2D, Point3D


class ModelPoint2D(Point2D):
    """
    Point in 2D, part of a 2D model
    @field _s - local coordinate of the point relative to the control grid's base S
    @field _t - local coordinate of the point relative to the control grid's base T
    @field _h - horizontal bilinear interpolant of the point relative to the cell it is in
    @field _v - vertical bilinear interpolant of the point relative to the cell it is in
    @field _cell_x - x coord of the left lower vertex of the cell that the point is in
    @field _cell_y - y coord of the left lower vertex of the cell that the point is in
    """
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
    """
    Point in 3D, part of a 3D model
    @field _s - local coordinate of the point relative to the control grid's base S
    @field _t - local coordinate of the point relative to the control grid's base T
    @field _u - local coordinate of the point relative to the control grid's base U
    """
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
