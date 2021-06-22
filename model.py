from point import Point2D, Point3D
from model_point import ModelPoint2D, ModelPoint3D
from control_grid import ControlGrid2D, ControlGrid3D


class Model2D:
    def __init__(self, *pts):
        self._vertices = []
        for point in pts:
            if not isinstance(point, Point2D):
                raise TypeError("Expected: %s, received: %s", Point2D.__name__, type(point).__name__)
            modelPoint = ModelPoint2D(point.x, point.y)
            self._vertices.append(modelPoint)
        self._min_x = min(list(map(lambda p: p.x, self._vertices)))
        self._min_y = min(list(map(lambda p: p.y, self._vertices)))
        self._max_x = max(list(map(lambda p: p.x, self._vertices)))
        self._max_y = max(list(map(lambda p: p.y, self._vertices)))

    @property
    def min_x(self):
        return self._min_x

    @property
    def min_y(self):
        return self._min_y

    @property
    def max_x(self):
        return self._max_x

    @property
    def max_y(self):
        return self._max_y

    @property
    def vertices(self):
        return self._vertices


class Model3D:
    def __init__(self, *pts):
        self._vertices = []
        for point in pts:
            if not isinstance(point, Point3D):
                raise TypeError("Expected: %s, received: %s", Point3D.__name__, type(point).__name__)
            modelPoint = ModelPoint3D(point.x, point.y, point.z)
            self._vertices.append(modelPoint)
        self._min_x = min(list(map(lambda p: p.x, self._vertices)))
        self._min_y = min(list(map(lambda p: p.y, self._vertices)))
        self._min_z = min(list(map(lambda p: p.z, self._vertices)))
        self._max_x = max(list(map(lambda p: p.x, self._vertices)))
        self._max_y = max(list(map(lambda p: p.y, self._vertices)))
        self._max_z = max(list(map(lambda p: p.z, self._vertices)))

    @property
    def min_x(self):
        return self._min_x

    @property
    def min_y(self):
        return self._min_y

    @property
    def min_z(self):
        return self._min_z

    @property
    def max_x(self):
        return self._max_x

    @property
    def max_y(self):
        return self._max_y

    @property
    def max_z(self):
        return self._max_z

    @property
    def vertices(self):
        return self._vertices
