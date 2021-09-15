import math

import numpy as np

import file_reader
from model_point import ModelPoint2D, ModelPoint3D
from point import Point2D, Point3D

# Define constants
MODEL_DIR = 'models/'
SPHERE_OBJ = 'sphere.obj'
TORUS_OBJ = 'torus.obj'


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

    @staticmethod
    def circle(radius):
        angle_values = np.linspace(0, 2 * math.pi, 100)
        points = []
        for angle in angle_values:
            x = radius * math.cos(angle)
            y = radius * math.sin(angle)
            points.append(Point2D(x, y))
        return Model2D(*points)

    @staticmethod
    def flower():
        angle_values = np.linspace(0, 2 * math.pi, 100)
        points = []
        a = 4
        k = 2
        for angle in angle_values:
            x = a * math.cos(k * angle) * math.cos(angle)
            y = a * math.cos(k * angle) * math.sin(angle)
            points.append(Point2D(x, y))
        return Model2D(*points)


class Model3D:
    def __init__(self, *pts):
        self._vertices = []
        self._triangles = []
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

    @property
    def triangles(self):
        return self._triangles

    @triangles.setter
    def triangles(self, t):
        self._triangles = t

    @staticmethod
    def sphere():
        vertices, triangles = file_reader.read_obj(MODEL_DIR + SPHERE_OBJ)
        points = []
        for vertex in vertices:
            points.append(Point3D(vertex[0], vertex[1], vertex[2]))
        model = Model3D(*points)
        model.triangles = triangles
        return model

    @staticmethod
    def torus():
        vertices, triangles = file_reader.read_obj(MODEL_DIR + TORUS_OBJ)
        points = []
        for vertex in vertices:
            points.append(Point3D(vertex[0], vertex[1], vertex[2]))
        model = Model3D(*points)
        model.triangles = triangles
        return model
