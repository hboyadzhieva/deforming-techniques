import math

import numpy as np
from mpl_toolkits.mplot3d import proj3d

import calculation as calc
from control_grid import ControlGrid2D, ControlGrid3D
from model import Model2D, Model3D
from point import Point2D, Point3D


# FIXME optimize comment style and keep consistent
def xs_ys_from_vertex_list(vertex_list):
    """ Return lists of coordinate values from a given vertex list

    Parameters:
    vertex_list: list
        list of Point2D instances

    Returns:
    list:
       list of x coordinates of the points
    list:
       list of y coordinates of the points
    """
    return list(map(lambda p: p.x, vertex_list)), list(map(lambda p: p.y, vertex_list))


def xs_ys_zs_from_vertex_list(vertex_list):
    """ Return lists of coordinate values from a given vertex list

    Parameters
    ----------
    vertex_list: list
        list of Point3D instances

    Returns
    --------
    list:
       list of x coordinates of the points
    list:
       list of y coordinates of the points
    list:
       list of y coordinates of the points
    """
    return np.array(list(map(lambda p: p.x, vertex_list))), np.array(list(map(lambda p: p.y, vertex_list))), np.array(
        list(map(lambda p: p.z, vertex_list)))


class GridDeformator2D:

    def __init__(self, model, control_points_x=4, control_points_y=4, offset_grid_model=2, offset_mouse_touch=0.2):
        if not isinstance(model, Model2D):
            raise TypeError("Type error. Expected: %s, received: %s", Model2D.__name__, type(model).__name__)
        self._model = model

        # Interactor settings
        self._offset_mouse_touch = offset_mouse_touch
        self._vertex_on_move = None

        # Choose center, S and T vectors so that the grid covers completely the model and leaves a margin of "offset"
        # units around the model. Initialize control grid with given number of control points per direction.
        center = Point2D(self._model.min_x - offset_grid_model, self._model.min_y - offset_grid_model)
        S = Point2D(self._model.max_x - center.x + offset_grid_model, 0)
        T = Point2D(0, self._model.max_y - center.y + offset_grid_model)
        self._grid = ControlGrid2D(control_points_x, control_points_y, center, S, T)

        self.setup_model_grid_params()

    def setup_model_grid_params(self):

        # Find parameters s,t,h and v of each vertex from the model according to the grid and grid's control points
        # positions
        for model_point in self._model.vertices:
            # calculate projections of Vertex P on S and T vectors
            pp0 = np.subtract(model_point.to_numpy_array(), self._grid.center.to_numpy_array())
            model_point.s = calc.scalar_projection(pp0, self._grid.S.to_numpy_array())
            model_point.t = calc.scalar_projection(pp0, self._grid.T.to_numpy_array())

            # FIXME move these calculations to Model, have fields for cell reference
            # calculate bilinear interpolants of the point relevant to the
            # 4 control points that define the cell in which the vertex is positioned
            si = (model_point.s * (self._grid.count_x_points - 1)) / np.linalg.norm(self._grid.S.to_numpy_array())
            tj = (model_point.t * (self._grid.count_y_points - 1)) / np.linalg.norm(self._grid.T.to_numpy_array())
            # largest integer smaller than si
            # for points on the lower or lefter boundary, consider them part of their upper/righer cells
            i = 0 if math.ceil(si) - 1 < 0 else math.ceil(si) - 1
            j = 0 if math.ceil(tj) - 1 < 0 else math.ceil(tj) - 1
            model_point.h = (model_point.x - self._grid.control_points[i][j].x) / (
                    self._grid.control_points[i + 1][j].x - self._grid.control_points[i][j].x)
            model_point.v = (model_point.y - self._grid.control_points[i][j].y) / (
                    self._grid.control_points[i][j + 1].y - self._grid.control_points[i][j].y)

    def plot_grid(self, ax):
        """
        Plot grid for 2D Grid Deformation
        Plot lines along x-axis with control points
        Plot lines along y-axis with control points
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
        """Plot model with vertices positions according to 2D Grid Deformation of the Grid

        Parameters
        ---------
        ax: Axes to plot on
        """
        for model_point in self._model.vertices:
            # FIXME use model point cell reference fields instead of calculating again
            si = (model_point.s * (self._grid.count_x_points - 1)) / np.linalg.norm(self._grid.S.to_numpy_array())
            tj = (model_point.t * (self._grid.count_y_points - 1)) / np.linalg.norm(self._grid.T.to_numpy_array())
            i = 0 if math.ceil(si) - 1 < 0 else math.ceil(si) - 1
            j = 0 if math.ceil(tj) - 1 < 0 else math.ceil(tj) - 1
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
            if calc.distance_2d(point_clicked, control_point) < self._offset_mouse_touch:
                self._vertex_on_move = control_point
                break

    def on_mouse_motion(self, event):
        if self._vertex_on_move is not None:
            x = event.xdata
            y = event.ydata
            self._vertex_on_move.x = x
            self._vertex_on_move.y = y

    def on_mouse_up(self, event):
        if self._vertex_on_move is not None:
            self._vertex_on_move = None


def coords_from_click_on_line(xd, yd, p0, p1, ax):
    # FIXME move this to appropriate place
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


class FreeFormDeformator():
    # FIXME add consistent comments on class and all methods
    def __init__(self, model, control_points_x=3, control_points_y=3, control_points_z=3, offset_grid_model=0,
                 offset_mouse_event=0.2):
        # FIXME think about consistency with public and private properties
        if not isinstance(model, Model3D):
            raise TypeError("Type error. Expected: %s, received: %s", Model3D.__name__, type(model).__name__)
        self._model = model

        # Setup Interactor parameters
        self.offset_mouse_event = offset_mouse_event
        self.offset_grid_model = offset_grid_model
        self._vertex_selected = None
        self._line_picked = None
        self._original_vertex_selected = None
        self.d = 0

        # Choose center, S,T,U vectors so that the grid covers completely the model and leaves a margin of "offset"
        # units around the model. Initialize control grid with given number of control points per direction.
        center = Point3D(self._model.min_x - offset_grid_model, self._model.min_y - offset_grid_model,
                         self._model.min_z - offset_grid_model)
        S = Point3D(self._model.max_x - center.x + offset_grid_model, 0, 0)
        T = Point3D(0, self._model.max_y - center.y + offset_grid_model, 0)
        U = Point3D(0, 0, self._model.max_z - center.z + offset_grid_model)
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

    def plot_grid(self, ax):
        """
        Plot grid for Free-Form Deformation
        Plot lines along x-axis with control points
        Plot lines along y-axis with control points
        Plot lines along z-axis with control points
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
        for model_point in self._model.vertices:
            new_vertex = calc.bezier_volume(self._grid, model_point.s, model_point.t, model_point.u)
            model_point.x = new_vertex[0]
            model_point.y = new_vertex[1]
            model_point.z = new_vertex[2]
        xs, ys, zs = xs_ys_zs_from_vertex_list(self._model.vertices)
        # FIXME use property for triangles
        ax.plot_trisurf(xs, ys, self._model._triangles, Z=zs, shade=True, color='white')

    def on_mouse_down(self, event, ax):
        # after moving the arrows, user clicks on another point from the canvas
        if self._vertex_selected is not None:
            self._vertex_selected = None
            self._original_vertex_selected = None
            self._line_picked = None
        if self._line_picked is not None and self._vertex_selected is None:
            xclick, yclick = event.xdata, event.ydata
            n = len(self._line_picked.get_data_3d()) - 1
            p0, p1 = np.transpose(self._line_picked.get_data_3d())[0], np.transpose(self._line_picked.get_data_3d())[n]
            pp0, pp1 = proj3d.proj_points((p0, p1), ax.M)
            xm, ym, zm = coords_from_click_on_line(xclick, yclick, pp0, pp1, ax)
            point_clicked = Point3D(xm, ym, zm)
            for control_point in self._grid.flat_control_points():
                if calc.distance_3d(point_clicked, control_point) < self.offset_mouse_event:
                    self._vertex_selected = control_point
                    self._original_vertex_selected = Point3D(control_point.x, control_point.y, control_point.z)
                    break

    def on_pick(self, event):
        if self._vertex_selected is None:
            self._line_picked = event.artist

    def add_arrows(self, ax):
        # FIXME add arrows as field of the class in init method
        self.x_arrow = ax.quiver(
            self._vertex_selected.x, self._vertex_selected.y, self._vertex_selected.z,  # <-- starting point of vector
            1, 0, 0,  # <-- directions of vector
            color='red', alpha=.8, lw=1, picker=True
        )
        self.y_arrow = ax.quiver(
            self._vertex_selected.x, self._vertex_selected.y, self._vertex_selected.z,
            # <-- starting point of vector
            0, 1, 0,  # <-- directions of vector
            color='blue', alpha=.8, lw=1, picker=True
        )
        self.z_arrow = ax.quiver(
            self._vertex_selected.x, self._vertex_selected.y, self._vertex_selected.z,
            # <-- starting point of vector
            0, 0, 1,  # <-- directions of vector
            color='green', alpha=.8, lw=1, picker=True
        )
