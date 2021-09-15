"""
Contains classes that reprsent the models that are subject to deformation
"""
import math

import numpy as np

import file_reader
import constants as cnt
from model_point import ModelPoint2D, ModelPoint3D
from point import Point2D, Point3D


class Model2D:
    """
    Model in 2D space
    @fields _min_x, _min_y: represent the minimum coordinates of vertices from the model
    @fields _max_x, _max_y: represent the maximum coordinates of vertices from the model
    @field _vertices: a list of all vertices (model_point.ModelPoint2D) that comprise the model
    @method circle: create a circle (Model2D)
    @method flower: create a flower (Model2D)
    """
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
        """
        Populate model with vertices on a circle.

        :param radius: the radius of the circle
        :return: circle model (Model2D)
        """
        angle_values = np.linspace(0, 2 * math.pi, 100)
        points = []
        for angle in angle_values:
            x = radius * math.cos(angle)
            y = radius * math.sin(angle)
            points.append(Point2D(x, y))
        return Model2D(*points)

    @staticmethod
    def flower():
        """
        Populate model with vertices on a flower.

        :return: flower model (Model2D)
        """
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
    """
    Model in 3D space
    @fields _min_x, _min_y, _min_z: represent the minimum coordinates of vertices from the model
    @fields _max_x, _max_y, _max_z: represent the maximum coordinates of vertices from the model
    @field _vertices: a list of all vertices (model_point.ModelPoint3D) that comprise the model
    @field _triangles: a list of all faces of the model, a face is a list of 3 verticess
    @method sphere: create a sphere (Model3D)
    @method torus: create a torus (Model3D)
    """
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
        """
        Populate model with vertices and triangles from file sphere.obj

        :return: sphere model (Model3D)
        """
        vertices, triangles = file_reader.read_obj(cnt.MODEL_DIR + cnt.SPHERE_OBJ)
        points = []
        for vertex in vertices:
            points.append(Point3D(vertex[0], vertex[1], vertex[2]))
        model = Model3D(*points)
        model.triangles = triangles
        return model

    @staticmethod
    def torus():
        """
        Populate model with vertices and triangles from file torus.obj

        :return: torus model (Model3D)
        """
        vertices, triangles = file_reader.read_obj(cnt.MODEL_DIR + cnt.TORUS_OBJ)
        points = []
        for vertex in vertices:
            points.append(Point3D(vertex[0], vertex[1], vertex[2]))
        model = Model3D(*points)
        model.triangles = triangles
        return model
