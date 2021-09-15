"""
A point in n-dim space
"""
import numpy


class Point2D:
    """
    Point in 2D space with x and y coordinatess
    """
    def __init__(self, x, y):
        self._x = x
        self._y = y

    @property
    def x(self):
        return self._x

    @x.setter
    def x(self, x):
        self._x = x

    @property
    def y(self):
        return self._y

    @y.setter
    def y(self, y):
        self._y = y

    def to_numpy_array(self):
        return numpy.array([self.x, self.y])


class Point3D:
    """
    Point in 3D space with x, y and z coordinatess
    """
    def __init__(self, x, y, z):
        self._x = x
        self._y = y
        self._z = z

    @property
    def x(self):
        return self._x

    @x.setter
    def x(self, x):
        self._x = x

    @property
    def y(self):
        return self._y

    @y.setter
    def y(self, y):
        self._y = y

    @property
    def z(self):
        return self._z

    @z.setter
    def z(self, z):
        self._z = z

    def to_numpy_array(self):
        return numpy.array([self.x, self.y, self.z])
