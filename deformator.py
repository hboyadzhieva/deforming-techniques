import math
from abc import ABC, abstractmethod

import numpy as np

import calculation as calc
from control_grid import ControlGrid2D, ControlGrid3D
from model import Model2D, Model3D
from point import Point2D, Point3D


class Deformator(ABC):
    def start(self):
        pass


class GridDeformator2D(Deformator):

    def __init__(self, model, control_points_x=4, control_points_y=4, offset=2):
        if not isinstance(model, Model2D):
            raise TypeError("Type error. Expected: %s, received: %s", Model2D.__name__, type(model).__name__)
        self._model = model

        # Choose center, S and T vectors so that the grid covers completely the model and leaves a margin of "offset"
        # units around the model. Initialize control grid with given number of control points per direction.
        center = Point2D(self._model.min_x - offset, self._model.min_y - offset)
        S = Point2D(self._model.max_x - center.x + offset, 0)
        T = Point2D(0, self._model.max_y - center.y + offset)
        self._grid = ControlGrid2D(control_points_x, control_points_y, center, S, T)

        self.setup_model_grid_params()

    def setup_model_grid_params(self):

        # Find parameters s,t,h and v of each vertex from the model according to the grid and grid's control points
        # positions
        for model_point in self._model.vertices:
            # calculate projections of Vertex P on S and T vectors
            pp0 = np.subtract(model_point.to_numpy_array(), self._grid.center.to_numpy_array())
            model_point.s = calc.vector_projection(pp0, self._grid.S.to_numpy_array())
            model_point.t = calc.vector_projection(pp0, self._grid.T.to_numpy_array())

            # calculate bilinear interpolants of the point relevant to the
            # 4 control points that define the cell in which the vertex is positioned
            right = math.ceil(model_point.s)
            up = math.ceil(model_point.t)
            model_point.h = (model_point.x - self._grid.control_points[right - 1][up - 1].x) / (
                    self._grid.control_points[right][up - 1].x - self._grid.control_points[right - 1][up - 1].x)
            model_point.v = (model_point.y - self._grid.control_points[right][up - 1].y) / (
                    self._grid.control_points[right][up].y - self._grid.control_points[right][up - 1].y)

    def start(self):

        print("Started 2D Grid Deformator")


class FreeFormDeformator(Deformator):
    def __init__(self, model, control_points_x=3, control_points_y=3, control_points_z=3, offset=2):
        if not isinstance(model, Model3D):
            raise TypeError("Type error. Expected: %s, received: %s", Model3D.__name__, type(model).__name__)
        self._model = model

        # Choose center, S,T,U vectors so that the grid covers completely the model and leaves a margin of "offset"
        # units around the model. Initialize control grid with given number of control points per direction.
        center = Point3D(self._model.min_x - offset, self._model.min_y - offset, self._model.min_z - offset)
        S = Point3D(self._model.max_x - center.x + offset, 0, 0)
        T = Point3D(0, self._model.max_y - center.y + offset, 0)
        U = Point3D(0, 0, self._model.max_z - center.z + offset)
        self._grid = ControlGrid3D(control_points_x, control_points_y, control_points_z, center, S, T, U)

        self.setup_model_grid_params()

    def setup_model_grid_params(self):
        for model_point in self._model.vertices:
            pp0 = np.subtract(model_point.to_numpy_array(), self._grid.center.to_numpy_array())
            model_point.s = calc.trilinear_interpolant(self._grid.T.to_numpy_array(), self._grid.U.to_numpy_array(),
                                                       pp0, self._grid.S.to_numpy_array())
            model_point.t = calc.trilinear_interpolant(self._grid.U.to_numpy_array(), self._grid.S.to_numpy_array(),
                                                       pp0, self._grid.T.to_numpy_array())
            model_point.u = calc.trilinear_interpolant(self._grid.S.to_numpy_array(), self._grid.T.to_numpy_array(),
                                                       pp0, self._grid.U.to_numpy_array())
