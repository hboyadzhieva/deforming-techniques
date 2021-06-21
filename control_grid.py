import numpy as np

from point import Point2D


class ControlGrid2D:
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

    @count_x_points.setter
    def count_x_points(self, count_x_points):
        self._count_x_points = count_x_points

    @property
    def count_y_points(self):
        return self._count_y_points

    @count_y_points.setter
    def count_y_points(self, count_y_points):
        self._count_y_points = count_y_points

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

    @control_points.setter
    def control_points(self, control_points):
        self._control_points = control_points

    def init_control_points(self):
        for i in range(0, self.count_x_points):
            buffer = []
            for j in range(0, self.count_y_points):
                # Pij = P0 + (i/n) * S + (j/m) * T,
                # where n and m are the number of segments obtained from the control points
                Pij = np.add(self._center.to_numpy_array(),
                             np.add((i / (self.count_x_points - 1)) * self._S.to_numpy_array(),
                                    (j / (self.count_y_points - 1)) * self._T.to_numpy_array()))
                buffer.append(Point2D(Pij[0], Pij[1]))
            self.control_points.append(buffer)

    def flat_control_points(self):
        result = []
        for i in self._control_points:
            for j in self._control_points[i]:
                result.append(j)
        return result


class ControlGrid3D:
    def __init__(self, count_x_points, count_y_points, count_z_points):
        self._count_x_points = count_x_points
        self._count_y_points = count_y_points
        self._count_z_points = count_z_points
        self._S = None
        self._T = None
        self._U = None
        self._center = None
        self.control_points = [[[]]]

    @property
    def count_x_points(self):
        return self._count_x_points

    @count_x_points.setter
    def count_x_points(self, count_x_points):
        self._count_x_points = count_x_points

    @property
    def count_y_points(self):
        return self._count_y_points

    @count_y_points.setter
    def count_y_points(self, count_y_points):
        self._count_y_points = count_y_points

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

    @control_points.setter
    def control_points(self, control_points):
        self._control_points = control_points

    def flat_control_points(self):
        result = []
        for i in self._control_points:
            for j in self._control_points[i]:
                for k in self._control_points[j]:
                    result.append(k)
        return result
