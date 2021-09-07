"""
Deformator is responsible for calculating the positions of the model vertices
related to the current grid's control points' positions.
"""

import math
import numpy as np
import calculation as calc
from mpl_toolkits.mplot3d import proj3d
from control_grid import ControlGrid2D, ControlGrid3D
from model import Model2D, Model3D
from point import Point2D, Point3D


class GridDeformator2D:
    """
    Simulator of 2D Grid Deformation.\n
    @const OFFSET_MOUSE_TOUCH: radius that defines area around mouse click, that counts as clicked on (float)\n
    @const OFFSET_GRID_MODEL: distance between model endpoints and control grid frame (float)\n
    @field vertex_selected: vertex that is currently being displaced by the user (point.Point2D)\n
    @method plot_grid, plot_model\n
    @method on_mouse_down, on_mouse_up, on_mouse_motion\n
    """

    OFFSET_MOUSE_TOUCH = 0.2
    OFFSET_GRID_MODEL = 2

    def __init__(self, model, control_points_x=4, control_points_y=4):
        if not isinstance(model, Model2D):
            raise TypeError("Type error. Expected: %s, received: %s", Model2D.__name__, type(model).__name__)
        self.vertex_selected = None
        self._model = model
        center = self._setup_grid_center()
        S, T = self._setup_grid_bases(center)
        self._grid = ControlGrid2D(control_points_x, control_points_y, center, S, T)
        self._setup_model_params()

    def plot_grid(self, ax):
        """
        Plot grid for 2D Grid Deformation.
        Plot lines along x-axis with control points.
        Plot lines along y-axis with control points.

        @:param ax: axes to plot on (matplotlib.Axes)
        """
        for j in range(0, self._grid.count_y_points):
            x_axis_line = []
            for i in range(0, self._grid.count_x_points):
                x_axis_line.append(self._grid.control_points[i][j])
            xs, ys = xs_ys_from_vertex_list(x_axis_line)
            ax.plot(xs, ys, 'k.-')
            x_axis_line.clear()

        for i in range(0, self._grid.count_x_points):
            y_axis_line = []
            for j in range(0, self._grid.count_y_points):
                y_axis_line.append(self._grid.control_points[i][j])
            xs, ys = xs_ys_from_vertex_list(y_axis_line)
            ax.plot(xs, ys, 'k.-')
            y_axis_line.clear()

    def plot_model(self, ax):
        """
        Plot model using bilinear interpolation to determine vertex position

        @:param ax: axes to plot on (matplotlib.Axes)
        """
        for model_point in self._model.vertices:
            i = model_point.cell_x
            j = model_point.cell_y
            new_vertex = calc.bilinear_interpolation(self._grid.control_points[i][j].to_numpy_array(),
                                                     self._grid.control_points[i + 1][j].to_numpy_array(),
                                                     self._grid.control_points[i + 1][j + 1].to_numpy_array(),
                                                     self._grid.control_points[i][j + 1].to_numpy_array(),
                                                     model_point.h, model_point.v)
            model_point.x = new_vertex[0]
            model_point.y = new_vertex[1]

        xs, ys = xs_ys_from_vertex_list(self._model.vertices)
        ax.plot(xs, ys, 'g-')

    def on_mouse_down(self, event):
        x = event.xdata
        y = event.ydata
        point_clicked = Point2D(x, y)
        for control_point in self._grid.flat_control_points():
            if calc.distance_2d(point_clicked, control_point) < self.OFFSET_MOUSE_TOUCH:
                self.vertex_selected = control_point
                break

    def on_mouse_motion(self, event):
        if self.vertex_selected is not None:
            x = event.xdata
            y = event.ydata
            self.vertex_selected.x = x
            self.vertex_selected.y = y

    def on_mouse_up(self, event):
        if self.vertex_selected is not None:
            self.vertex_selected = None

    def _setup_grid_center(self):
        """
        Choose appropriate grid center position so that all model vertices are inside the grid

        @:returns vertex, center of control grid (point.Point2D)
        """
        center = Point2D(self._model.min_x - self.OFFSET_GRID_MODEL, self._model.min_y - self.OFFSET_GRID_MODEL)
        return center

    def _setup_grid_bases(self, center):
        """
        Choose appropriate base vector positions so that all model vertices are inside the grid

        @:returns vectors representing the bases of the control grid (point.Point2D, point.Point2D)
        """
        S = Point2D(self._model.max_x - center.x + self.OFFSET_GRID_MODEL, 0)
        T = Point2D(0, self._model.max_y - center.y + self.OFFSET_GRID_MODEL)
        return S, T

    def _setup_model_params(self):
        """
        Setup parameters s,t,h,v,cell_x and cell_y of each vertex from
        the model according to the grid and grid's control points positions.
        """
        for model_point in self._model.vertices:
            pp0 = np.subtract(model_point.to_numpy_array(), self._grid.center.to_numpy_array())
            model_point.s = calc.scalar_projection(pp0, self._grid.S.to_numpy_array())
            model_point.t = calc.scalar_projection(pp0, self._grid.T.to_numpy_array())

            si = (model_point.s * (self._grid.count_x_points - 1)) / np.linalg.norm(self._grid.S.to_numpy_array())
            tj = (model_point.t * (self._grid.count_y_points - 1)) / np.linalg.norm(self._grid.T.to_numpy_array())
            # largest integer smaller than si
            # for points on the lower or lefter boundary, consider them part of their upper/righer cells
            i = model_point.cell_x = 0 if math.ceil(si) - 1 < 0 else math.ceil(si) - 1
            j = model_point.cell_y = 0 if math.ceil(tj) - 1 < 0 else math.ceil(tj) - 1
            model_point.h = (model_point.x - self._grid.control_points[i][j].x) / (
                    self._grid.control_points[i + 1][j].x - self._grid.control_points[i][j].x)
            model_point.v = (model_point.y - self._grid.control_points[i][j].y) / (
                    self._grid.control_points[i][j + 1].y - self._grid.control_points[i][j].y)


class FreeFormDeformator:
    """
    Simulator of 3D Grid Deformation.\n
    @const OFFSET_MOUSE_TOUCH: radius that defines area around mouse click, that counts as clicked on (float)\n
    @const OFFSET_GRID_MODEL: distance between model endpoints and control grid frame (float)\n
    @field vertex_selected: vertex that is currently being displaced by the user (point.Point3D)\n
    @field original_vertex_selected: original position of vertex_selected before being displaced (point.Point3D)\n
    @method plot_grid, plot_model\n
    @method on_mouse_down, on_mouse_up, on_mouse_motion\n
    """

    OFFSET_MOUSE_TOUCH = 0.2
    OFFSET_GRID_MODEL = 0.1

    def __init__(self, model, control_points_x=3, control_points_y=3, control_points_z=3):
        if not isinstance(model, Model3D):
            raise TypeError("Type error. Expected: %s, received: %s", Model3D.__name__, type(model).__name__)
        self._model = model
        self.vertex_selected = None
        self.original_vertex_selected = None
        self._line_picked = None
        self.d = 0
        center = self._setup_grid_center()
        S, T, U = self._setup_grid_bases(center)
        self._grid = ControlGrid3D(control_points_x, control_points_y, control_points_z, center, S, T, U)
        self._setup_model_grid_params()

    def plot_grid(self, ax):
        """
        Plot grid for Free-Form Deformation \n
        Plot lines along x-axis with control points \n
        Plot lines along y-axis with control points \n
        Plot lines along z-axis with control points \n

        @:param ax: axes to plot on (matplotlib.Axes3D)
        """

        for j in range(0, self._grid.count_y_points):
            for k in range(0, self._grid.count_z_points):
                x_axis_line = []
                for i in range(0, self._grid.count_x_points):
                    x_axis_line.append(self._grid.control_points[i][j][k])
                xs, ys, zs = xs_ys_zs_from_vertex_list(x_axis_line)
                ax.plot(xs, ys, zs, 'k.-', picker=5)
                x_axis_line.clear()

        for i in range(0, self._grid.count_x_points):
            for k in range(0, self._grid.count_z_points):
                y_axis_line = []
                for j in range(0, self._grid.count_y_points):
                    y_axis_line.append(self._grid.control_points[i][j][k])
                xs, ys, zs = xs_ys_zs_from_vertex_list(y_axis_line)
                ax.plot(xs, ys, zs, 'k.-')
                y_axis_line.clear()

        for i in range(0, self._grid.count_x_points):
            for j in range(0, self._grid.count_y_points):
                z_axis_line = []
                for k in range(0, self._grid.count_z_points):
                    z_axis_line.append((self._grid.control_points[i][j][k]))
                xs, ys, zs = xs_ys_zs_from_vertex_list(z_axis_line)
                ax.plot(xs, ys, zs, 'k.-')
                z_axis_line.clear()

    def plot_model(self, ax):
        """
        Plot model using bezier interpolation to determine vertex position

        @:param ax: axes to plot on (matplotlib.Axes3D)
        """
        for model_point in self._model.vertices:
            new_vertex = calc.point_position_in_bezier_volume(self._grid, model_point.s, model_point.t, model_point.u)
            model_point.x = new_vertex[0]
            model_point.y = new_vertex[1]
            model_point.z = new_vertex[2]
        xs, ys, zs = xs_ys_zs_from_vertex_list(self._model.vertices)
        ax.plot_trisurf(xs, ys, self._model.triangles, Z=zs, shade=True, color='white')

    def on_mouse_down(self, event, ax):
        if self.vertex_selected is not None:
            self.vertex_selected = None
            self.original_vertex_selected = None
            self._line_picked = None
        if self._line_picked is not None and self.vertex_selected is None:
            xclick, yclick = event.xdata, event.ydata
            n = len(self._line_picked.get_data_3d()) - 1
            p0, p1 = np.transpose(self._line_picked.get_data_3d())[0], np.transpose(self._line_picked.get_data_3d())[n]
            pp0, pp1 = proj3d.proj_points((p0, p1), ax.M)
            xm, ym, zm = coords_from_click_on_line(xclick, yclick, pp0, pp1, ax)
            point_clicked = Point3D(xm, ym, zm)
            for control_point in self._grid.flat_control_points():
                if calc.distance_3d(point_clicked, control_point) < self.OFFSET_MOUSE_TOUCH:
                    self.vertex_selected = control_point
                    self.original_vertex_selected = Point3D(control_point.x, control_point.y, control_point.z)
                    break

    def on_pick(self, event):
        if self.vertex_selected is None:
            self._line_picked = event.artist

    def add_arrows(self, ax):
        """
        Add arrows around selected vertex to show coordinate directions.\n
        red arrow: x direction, blue arrow: y direction, green arrow: z direction
        """
        ax.quiver(
            self.vertex_selected.x, self.vertex_selected.y, self.vertex_selected.z,
            1, 0, 0, color='red', alpha=.8, lw=1, picker=True
        )
        ax.quiver(
            self.vertex_selected.x, self.vertex_selected.y, self.vertex_selected.z,
            0, 1, 0, color='blue', alpha=.8, lw=1, picker=True
        )
        ax.quiver(
            self.vertex_selected.x, self.vertex_selected.y, self.vertex_selected.z,
            0, 0, 1, color='green', alpha=.8, lw=1, picker=True
        )

    def _setup_grid_center(self):
        """
        Choose appropriate grid center position so that all model vertices are inside the grid

        @:returns vertex, center of control grid (point.Point3D)
        """
        center = Point3D(self._model.min_x - self.OFFSET_GRID_MODEL,
                         self._model.min_y - self.OFFSET_GRID_MODEL,
                         self._model.min_z - self.OFFSET_GRID_MODEL)
        return center

    def _setup_grid_bases(self, center):
        """
        Choose appropriate base vector positions so that all model vertices are inside the grid

        @:returns vectors representing the bases of the control grid (point.Point3D, point.Point3D, point.Point3D)
        """
        S = Point3D(self._model.max_x - center.x + self.OFFSET_GRID_MODEL, 0, 0)
        T = Point3D(0, self._model.max_y - center.y + self.OFFSET_GRID_MODEL, 0)
        U = Point3D(0, 0, self._model.max_z - center.z + self.OFFSET_GRID_MODEL)
        return S, T, U

    def _setup_model_grid_params(self):
        """
        Setup parameters s,t,u of each vertex from the model
        according to the grid and grid's control points positions.
        """
        for model_point in self._model.vertices:
            pp0 = np.subtract(model_point.to_numpy_array(), self._grid.center.to_numpy_array())
            model_point.s = calc.trilinear_interpolant(self._grid.T.to_numpy_array(), self._grid.U.to_numpy_array(),
                                                       pp0, self._grid.S.to_numpy_array())
            model_point.t = calc.trilinear_interpolant(self._grid.U.to_numpy_array(), self._grid.S.to_numpy_array(),
                                                       pp0, self._grid.T.to_numpy_array())
            model_point.u = calc.trilinear_interpolant(self._grid.S.to_numpy_array(), self._grid.T.to_numpy_array(),
                                                       pp0, self._grid.U.to_numpy_array())


def xs_ys_from_vertex_list(vertex_list):
    """
    Return lists of coordinate values from a given vertex list
    @:param vertex_list (list of point.Point2D)
    @:returns xs - x coordinates (list of floats), ys - y coordinates (list of floats)
    """
    return list(map(lambda p: p.x, vertex_list)), list(map(lambda p: p.y, vertex_list))


def xs_ys_zs_from_vertex_list(vertex_list):
    """
    Return lists of coordinate values from a given vertex list
    @:param vertex_list (list of point.Point2D)
    @:returns xs - x coords (list of floats), ys - y coords (list of floats), zs - z coords (list of floats)
    """
    return np.array(list(map(lambda p: p.x, vertex_list))), np.array(list(map(lambda p: p.y, vertex_list))), np.array(
        list(map(lambda p: p.z, vertex_list)))


def coords_from_click_on_line(xd, yd, p0, p1, ax):
    """
    Given the 2D view coordinates attempt to guess a 3D coordinate.
    Use method from axes3d but pass specific vertices for the edge.
    Finds a 3D coordinate on the edge with vertices p0 and p1
    """

    if ax.M is None:
        return ''

    # scale the z value to match
    x0, y0, z0 = p0
    x1, y1, z1 = p1
    d0 = np.hypot(x0 - xd, y0 - yd)
    d1 = np.hypot(x1 - xd, y1 - yd)
    dt = d0 + d1
    z = d1 / dt * z0 + d0 / dt * z1

    x, y, z = proj3d.inv_transform(xd, yd, z, ax.M)

    return x, y, z
