import re
from abc import ABC, abstractmethod

import matplotlib.pyplot as plt
import numpy as np
from mpl_toolkits import mplot3d
from mpl_toolkits.mplot3d import proj3d

import calculation as calc
from control_grid import ControlGrid2D, ControlGrid3D
from model import Model2D, Model3D
from point import Point2D, Point3D


class Deformator(ABC):
    def start(self):
        pass

    def register_mouse_events(self, fig):
        fig.canvas.mpl_connect('button_press_event', self.on_mouse_down)
        fig.canvas.mpl_connect('button_release_event', self.on_mouse_up)
        fig.canvas.mpl_connect('motion_notify_event', self.on_mouse_motion)
        fig.canvas.mpl_connect('pick_event', self.on_pick)

    def on_mouse_down(self, event):
        pass

    def on_mouse_up(self, event):
        pass

    def on_mouse_motion(self, event):
        pass

    def on_pick(self, event):
        pass

    def plot_grid(self):
        pass

    def plot_model(self):
        pass

    def plot(self):
        plt.cla()
        self.plot_grid()
        self.plot_model()
        plt.xlabel("x")
        plt.ylabel("y")
        plt.show()


def xs_ys_from_vertex_list(vertex_list):
    """ Return lists of coordinate values from a given vertex list

    Parameters
    ----------
    vertex_list: list
        list of Point2D instances

    Returns
    --------
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


class GridDeformator2D(Deformator):

    def __init__(self, model, control_points_x=4, control_points_y=4, offset_grid_model=2, offset_mouse_touch=0.2):
        if not isinstance(model, Model2D):
            raise TypeError("Type error. Expected: %s, received: %s", Model2D.__name__, type(model).__name__)
        self._model = model

        # Interactor settings
        self.ax = None
        self.fig = None
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
            model_point.s = calc.vector_projection(pp0, self._grid.S.to_numpy_array())
            model_point.t = calc.vector_projection(pp0, self._grid.T.to_numpy_array())

            # calculate bilinear interpolants of the point relevant to the
            # 4 control points that define the cell in which the vertex is positioned
            l = model_point.left()
            r = model_point.right()
            t = model_point.top()
            b = model_point.bottom()
            model_point.h = (model_point.x - self._grid.control_points[l][b].x) / (
                    self._grid.control_points[r][b].x - self._grid.control_points[l][b].x)
            model_point.v = (model_point.y - self._grid.control_points[r][b].y) / (
                    self._grid.control_points[r][t].y - self._grid.control_points[r][b].y)

    def start(self):
        print("Started 2D Grid Deformator")
        self.fig, self.ax = plt.subplots()
        self.register_mouse_events(self.fig)
        self.plot()

    def plot_grid(self):
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
            self.ax.plot(xs, ys, 'k.-')
            x_axis_line.clear()

        for i in range(0, self._grid.count_x_points):
            y_axis_line = []
            for j in range(0, self._grid.count_y_points):
                y_axis_line.append(self._grid.control_points[i][j])
            xs, ys = xs_ys_from_vertex_list(y_axis_line)
            self.ax.plot(xs, ys, 'k.-')
            y_axis_line.clear()

    def plot_model(self):
        """Plot model with vertices positions according to 2D Grid Deformation of the Grid

        Parameters
        ---------
        ax: Axes to plot on
        """
        for model_point in self._model.vertices:
            l = model_point.left()
            r = model_point.right()
            t = model_point.top()
            b = model_point.bottom()
            new_vertex = calc.bilinear_interpolation(self._grid.control_points[l][b].to_numpy_array(),
                                                     self._grid.control_points[r][b].to_numpy_array(),
                                                     self._grid.control_points[r][t].to_numpy_array(),
                                                     self._grid.control_points[l][t].to_numpy_array(),
                                                     model_point.h, model_point.v)
            model_point.x = new_vertex[0]
            model_point.y = new_vertex[1]

        xs, ys = xs_ys_from_vertex_list(self._model.vertices)
        self.ax.plot(xs, ys, 'go')

    def on_mouse_down(self, event):
        x = event.xdata
        y = event.ydata
        point_clicked = Point2D(x, y)
        for control_point in self._grid.flat_control_points():
            if calc.distance_2d(point_clicked, control_point) < self._offset_mouse_touch:
                print("clicked")
                self._vertex_on_move = control_point
                break

    def on_mouse_motion(self, event):
        if self._vertex_on_move is not None:
            x = event.xdata
            y = event.ydata
            self._vertex_on_move.x = x
            self._vertex_on_move.y = y
            self.plot()

    def on_mouse_up(self, event):
        if self._vertex_on_move is not None:
            self._vertex_on_move = None


def format_coord(xd, yd, p0, p1, ax):
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


class FreeFormDeformator(Deformator):
    def __init__(self, model, control_points_x=3, control_points_y=3, control_points_z=3, offset_grid_model=2,
                 offset_mouse_event=0.2):
        if not isinstance(model, Model3D):
            raise TypeError("Type error. Expected: %s, received: %s", Model3D.__name__, type(model).__name__)
        self._model = model

        # Setup Interactor parameters
        self.ax = None
        self.offset_mouse_event = offset_mouse_event
        self._vertex_on_move = None
        self._line_picked = None

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

    def start(self):
        print("Started Free-Form Deformator")
        self.fig = plt.figure()
        self.ax = plt.subplot(projection="3d")
        self.register_mouse_events(self.fig)
        self.plot()

    def plot_grid(self):
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
                self.ax.plot(xs, ys, zs, 'k.-', picker=5)
                x_axis_line.clear()

        for i in range(0, self._grid.count_x_points):
            for k in range(0, self._grid.count_z_points):
                y_axis_line = []
                for j in range(0, self._grid.count_y_points):
                    y_axis_line.append(self._grid.control_points[i][j][k])
                xs, ys, zs = xs_ys_zs_from_vertex_list(y_axis_line)
                self.ax.plot(xs, ys, zs, 'k.-')
                y_axis_line.clear()

        for i in range(0, self._grid.count_x_points):
            for j in range(0, self._grid.count_y_points):
                z_axis_line = []
                for k in range(0, self._grid.count_z_points):
                    z_axis_line.append((self._grid.control_points[i][j][k]))
                xs, ys, zs = xs_ys_zs_from_vertex_list(z_axis_line)
                self.ax.plot(xs, ys, zs, 'k.-')
                z_axis_line.clear()

    def plot_model(self):
        for model_point in self._model.vertices:
            new_vertex = calc.bezier_volume(self._grid, model_point.s, model_point.t, model_point.u)
            model_point.x = new_vertex[0]
            model_point.y = new_vertex[1]
            model_point.z = new_vertex[2]
        xs, ys, zs = xs_ys_zs_from_vertex_list(self._model.vertices)
        self.ax.plot(xs, ys, zs, 'go', picker=5)

    def on_mouse_down(self, event):
        if self._line_picked is not None:
            xclick, yclick = event.xdata, event.ydata
            pick_xs = self._line_picked.get_data_3d()[0]
            pick_ys = self._line_picked.get_data_3d()[1]
            pick_zs = self._line_picked.get_data_3d()[2]
            p0 = np.array([pick_xs[0], pick_ys[0], pick_zs[0]])
            p1 = np.array([pick_xs[1], pick_ys[1], pick_zs[1]])
            pp0, pp1 = proj3d.proj_points((p0, p1), self.ax.M)
            xm, ym, zm = format_coord(xclick, yclick, pp0, pp1, self.ax)
            print("mouse on line click (%s, %s, %s)" % (xm, ym, zm))
            str = self.ax.format_coord(xclick, yclick)

    def on_mouse_motion(self, event):
        pass

    def on_mouse_up(self, event):
        if self._vertex_on_move is not None:
            self._vertex_on_move = None

    def on_pick(self, event):
        print("picked")
        self._line_picked = event.artist
        data = self._line_picked.get_data()
        data3d = self._line_picked.get_data_3d()
        p1, p2, p3 = data3d
