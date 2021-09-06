"""
Control grids that contain the objects to be deformed.
Control grids are composed of control points.
Control grids are defined on a coordinate system with center and base vectors.
"""
import numpy as np
from point import Point2D, Point3D


class ControlGrid2D:
    """
    Control grid in 2D space.
    @field _count_x_points: number of control points in x direction (int)
    @field _count_y_points: number of control points in y direction (int)
    @field _S: base vector in x direction (point.Point2D), _T: base vector in y direction (point.Point2D)
    @field _center: position vector representing the center of the coordinate system (point.Point2D)
    @field _control_points: all control points composing the grid (vector[][] of point.Point2D)
    @method init_control_points: initialize _control_points and their positions
    @method flat_control_points: return flat 1D vector of the control points
    """

    def __init__(self, count_x_points, count_y_points, center, S, T):
        self._count_x_points = count_x_points
        self._count_y_points = count_y_points
        self._S = S
        self._T = T
        self._center = center
        self._control_points = []
        self.init_control_points()

    @property
    def count_x_points(self):
        return self._count_x_points

    @property
    def count_y_points(self):
        return self._count_y_points

    @property
    def S(self):
        return self._S

    @S.setter
    def S(self, S):
        self._S = S

    @property
    def T(self):
        return self._T

    @T.setter
    def T(self, T):
        self._T = T

    @property
    def center(self):
        return self._center

    @center.setter
    def center(self, center):
        self._center = center

    @property
    def control_points(self):
        return self._control_points

    def init_control_points(self):
        """
        Given number of control points in each coordinate direction
        and the base vectors and center of the coordinate system,
        find the global position of each control point
        and store the points in 2D vector _control_points.
        """
        for i in range(0, self.count_x_points):
            buffer = []
            for j in range(0, self.count_y_points):
                Pij = np.add(self._center.to_numpy_array(),
                             np.add((i / (self.count_x_points - 1)) * self.S.to_numpy_array(),
                                    (j / (self.count_y_points - 1)) * self.T.to_numpy_array()))
                buffer.append(Point2D(Pij[0], Pij[1]))
            self.control_points.append(buffer)

    def flat_control_points(self):
        """Return flat one dimensional array of all control points in the grid"""
        result = []
        for i in self._control_points:
            for j in i:
                result.append(j)
        return result


class ControlGrid3D:
    """
       Control grid in 2D space.
       @field _count_x_points: number of control points in x direction (int)
       @field _count_y_points: number of control points in y direction (int)
       @field _count_z_points: number of control points in z direction (int)
       @field _S: base vector in x direction (point.Point3D), _T: base vector in y direction (point.Point3D)
       @field _U: base vector in z direction (point.Point3D)
       @field _center: position vector representing the center of the coordinate system (point.Point3D)
       @field _control_points: all control points composing the grid (vector[][][] of point.point3D)
       @method init_control_points: initialize _control_points and their positions
       @method flat_control_points: return flat 1D vector of the control points
       """
    def __init__(self, count_x_points, count_y_points, count_z_points, center, S, T, U):
        self._count_x_points = count_x_points
        self._count_y_points = count_y_points
        self._count_z_points = count_z_points
        self._S = S
        self._T = T
        self._U = U
        self._center = center
        self._control_points = []
        self.init_control_points()

    @property
    def count_x_points(self):
        return self._count_x_points

    @property
    def count_y_points(self):
        return self._count_y_points

    @property
    def count_z_points(self):
        return self._count_z_points

    @property
    def S(self):
        return self._S

    @S.setter
    def S(self, S):
        self._S = S

    @property
    def T(self):
        return self._T

    @T.setter
    def T(self, T):
        self._T = T

    @property
    def U(self):
        return self._U

    @U.setter
    def U(self, U):
        self._U = U

    @property
    def center(self):
        return self._center

    @center.setter
    def center(self, center):
        self._center = center

    @property
    def control_points(self):
        return self._control_points

    def init_control_points(self):
        """
        Given number of control points in each coordinate direction
        and the base vectors and center of the coordinate system,
        find the global position of each control point
        and store the points in 3D vector _control_points.
        """
        for i in range(0, self.count_x_points):
            row_buffer = []
            for j in range(0, self.count_y_points):
                col_buffer = []
                for k in range(0, self.count_z_points):
                    Pij = np.add(self._center.to_numpy_array(),
                                 np.add((i / (self.count_x_points - 1)) * self.S.to_numpy_array(),
                                        np.add((j / (self.count_y_points - 1)) * self.T.to_numpy_array(),
                                               (k / (self.count_z_points - 1)) * self.U.to_numpy_array())))
                    col_buffer.append(Point3D(Pij[0], Pij[1], Pij[2]))
                row_buffer.append(col_buffer)
            self.control_points.append(row_buffer)

    def flat_control_points(self):
        """Return flat one dimensional array of all control points in the grid"""
        result = []
        for i in self._control_points:
            for j in i:
                for k in j:
                    result.append(k)
        return result
