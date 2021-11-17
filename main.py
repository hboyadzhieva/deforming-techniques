"""
Program to demonstrate the 2D Grid Deformation and Free Form Deformation
techiques for deforming objects in 2D and 3D space.
"""
import tkinter as tk
import constants as cnt
import matplotlib.pyplot as plt
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg

from deformator import GridDeformator2D, FreeFormDeformator
from model import Model2D, Model3D


class MainApplication:
    """
    Create Tkinter window to display the techniques simulations.
    Switch between the two deformator, different models and interact
    with the canvas.
    """
    def __init__(self, r):
        self.root = r

        self.figure_2d = plt.Figure(figsize=(6, 5), dpi=100)
        self.ax_2d = self.figure_2d.add_subplot(111)
        self.gd2d = None
        self.figure_3d = plt.Figure(figsize=(6, 5), dpi=100)
        self.ax_3d = self.figure_3d.add_subplot(111, projection='3d')
        self.ffd = None
        self.canvas = FigureCanvasTkAgg(self.figure_2d, self.root)

        self.radio_def_value = tk.IntVar()
        self.radio_gd2d_value = tk.IntVar()
        self.radio_ffd_value = tk.IntVar()
        self.dx, self.dy, self.dz = tk.DoubleVar(), tk.DoubleVar(), tk.DoubleVar()
        self.check_show_grid_value = tk.IntVar(value=1)

        self.figure_2d.canvas.mpl_connect(cnt.BTN_PRESS_EVENT, self.on_mouse_down_2d)
        self.figure_2d.canvas.mpl_connect(cnt.BTN_RELEASE_EVENT, self.on_mouse_up_2d)
        self.figure_2d.canvas.mpl_connect(cnt.BTN_MOTION_EVENT, self.on_mouse_motion_2d)
        self.figure_3d.canvas.mpl_connect(cnt.BTN_PRESS_EVENT, self.on_mouse_down_3d)
        self.figure_3d.canvas.mpl_connect(cnt.PICK_EVENT, self.on_pick_3d)

        self.radio_gd2d = tk.Radiobutton(root)
        self.radio_ffd = tk.Radiobutton(root)
        self.radio_gd2d_circle = tk.Radiobutton(root)
        self.radio_gd2d_flower = tk.Radiobutton(root)
        self.radio_ffd_sphere = tk.Radiobutton(root)
        self.radio_ffd_torus = tk.Radiobutton(root)
        self.x_slide = tk.Scale(root)
        self.y_slide = tk.Scale(root)
        self.z_slide = tk.Scale(root)
        self.check_show_grid = tk.Checkbutton(root)

        self.widgets_2d = [self.radio_gd2d_circle, self.radio_gd2d_flower]
        self.widgets_3d = [self.radio_ffd_sphere, self.radio_ffd_torus, self.x_slide, self.y_slide, self.z_slide]

    def configure_gui_elements(self):
        """
        Configure tk widgets
        """
        self.root.geometry('%dx%d+%d+%d' % (1000, 800, 0, 0))
        self.root.title(cnt.TITLE)
        self.radio_gd2d.configure(text=cnt.GD2D, variable=self.radio_def_value, value=1,
                                  command=lambda: self.choose_def(cnt.GD2D))
        self.radio_ffd.configure(text=cnt.FFD, variable=self.radio_def_value, value=2,
                                 command=lambda: self.choose_def(cnt.FFD))
        self.radio_gd2d_circle.configure(text=cnt.CIRCLE, variable=self.radio_gd2d_value, value=1,
                                         command=lambda: self.display_2d(cnt.CIRCLE))
        self.radio_gd2d_flower.configure(text=cnt.FLOWER, variable=self.radio_gd2d_value, value=2,
                                         command=lambda: self.display_2d(cnt.FLOWER))
        self.radio_ffd_sphere.configure(text=cnt.SPHERE, variable=self.radio_ffd_value, value=1,
                                        command=lambda: self.display_3d(cnt.SPHERE))
        self.radio_ffd_torus.configure(text=cnt.TORUS, variable=self.radio_ffd_value, value=2,
                                       command=lambda: self.display_3d(cnt.TORUS))
        self.x_slide.configure(from_=-2, to=2, resolution=0.01, orient=tk.HORIZONTAL, variable=self.dx,
                               command=lambda x: self.change_pos(cnt.X_AXIS), fg=cnt.COLOR_RED, label=cnt.X_AXIS,
                               state=cnt.STATE_DISABLED)
        self.y_slide.configure(from_=-2, to=2, resolution=0.01, orient=tk.HORIZONTAL, variable=self.dy,
                               command=lambda x: self.change_pos(cnt.Y_AXIS), fg=cnt.COLOR_BLUE, label=cnt.Y_AXIS,
                               state=cnt.STATE_DISABLED)
        self.z_slide.configure(from_=-2, to=2, resolution=0.01, orient=tk.HORIZONTAL, variable=self.dz,
                               command=lambda x: self.change_pos(cnt.Z_AXIS), fg=cnt.COLOR_GREEN, label=cnt.Z_AXIS,
                               state=cnt.STATE_DISABLED)
        self.check_show_grid.configure(text=cnt.SHOW_GRID_TITLE, variable=self.check_show_grid_value, onvalue=1,
                                       offvalue=0,
                                       command=self.toggle_grid)

    def start(self):
        self.configure_gui_elements()
        self.radio_gd2d.pack()
        self.radio_ffd.pack()
        self.check_show_grid.pack()
        self.root.mainloop()

    def choose_def(self, deformator):
        """
        Switch between 2D Grid Deformator and Free Form Deformator simulators
        """
        self.canvas.get_tk_widget().pack_forget()
        if deformator == cnt.GD2D:
            for w_3d in self.widgets_3d:
                w_3d.pack_forget()
            for w_2d in self.widgets_2d:
                w_2d.pack()
        elif deformator == cnt.FFD:
            for w_2d in self.widgets_2d:
                w_2d.pack_forget()
            for w_3d in self.widgets_3d:
                w_3d.pack()

    def display_2d(self, model_choice):
        """
        Display 2D Grid Deformator with model of choice
        """
        self.canvas.get_tk_widget().pack_forget()
        self.ax_2d.cla()
        self.canvas = FigureCanvasTkAgg(self.figure_2d, root)
        self.canvas.get_tk_widget().pack()
        model = None
        if model_choice == cnt.CIRCLE:
            model = Model2D.circle(2)
        elif model_choice == cnt.FLOWER:
            model = Model2D.flower()
        self.gd2d = GridDeformator2D(model)
        self.gd2d.plot_model(self.ax_2d)
        if self.check_show_grid_value.get() == 1:
            self.gd2d.plot_grid(self.ax_2d)
        plt.show()

    def display_3d(self, model_choice):
        """
        Display Free Form Deformator with model of choice
        """
        self.canvas.get_tk_widget().pack_forget()
        self.ax_3d.cla()
        self.canvas = FigureCanvasTkAgg(self.figure_3d, root)
        self.canvas.get_tk_widget().pack()
        model = None
        if model_choice == cnt.SPHERE:
            model = Model3D.sphere()
        elif model_choice == cnt.TORUS:
            model = Model3D.torus()
        self.ffd = FreeFormDeformator(model)
        self.ffd.plot_model(self.ax_3d)
        self.ffd.plot_grid(self.ax_3d)
        plt.show()

    def change_pos(self, axis):
        """
        Displace the selected vertex along the axis together with the value on the slider
        """
        self.ax_3d.cla()
        vertex = self.ffd.vertex_selected
        original_vertex = self.ffd.original_vertex_selected
        if vertex is not None and original_vertex is not None:
            if axis == cnt.X_AXIS:
                vertex.x = original_vertex.x + self.dx.get()
            elif axis == cnt.Y_AXIS:
                vertex.y = original_vertex.y + self.dy.get()
            elif axis == cnt.Z_AXIS:
                vertex.z = original_vertex.z + self.dz.get()
            self.ffd.add_arrows(self.ax_3d)
        self.ffd.plot_model(self.ax_3d)
        if self.check_show_grid_value.get() == 1:
            self.ffd.plot_grid(self.ax_3d)
        self.figure_3d.canvas.draw()
        plt.show()

    def toggle_grid(self):
        """
        Show/hide the control grid
        """
        if self.radio_def_value.get() == 1:
            self.ax_2d.cla()
            self.gd2d.plot_model(self.ax_2d)
            if self.check_show_grid_value.get() == 1:
                self.gd2d.plot_grid(self.ax_2d)
            self.figure_2d.canvas.draw()
            plt.show()
        elif self.radio_def_value.get() == 2:
            self.ax_3d.cla()
            self.ffd.plot_model(self.ax_3d)
            if self.check_show_grid_value.get() == 1:
                self.ffd.plot_grid(self.ax_3d)
            self.figure_3d.canvas.draw()
            plt.show()

    def on_mouse_down_2d(self, event):
        self.gd2d.on_mouse_down(event)

    def on_mouse_up_2d(self, event):
        self.gd2d.on_mouse_up(event)

    def on_mouse_motion_2d(self, event):
        self.gd2d.on_mouse_motion(event)
        self.ax_2d.cla()
        self.gd2d.plot_model(self.ax_2d)
        if self.check_show_grid_value.get() == 1:
            self.gd2d.plot_grid(self.ax_2d)
        self.figure_2d.canvas.draw()
        plt.show()

    def on_mouse_down_3d(self, event):
        self.ffd.on_mouse_down(event, self.ax_3d)
        self.reset_sliders()
        if self.ffd.vertex_selected is not None:
            self.ax_3d.cla()
            self.ffd.add_arrows(self.ax_3d)
            self.ffd.plot_model(self.ax_3d)
            if self.check_show_grid_value.get() == 1:
                self.ffd.plot_grid(self.ax_3d)
            self.enable_sliders()
        else:
            self.ax_3d.cla()
            self.ffd.plot_model(self.ax_3d)
            if self.check_show_grid_value.get() == 1:
                self.ffd.plot_grid(self.ax_3d)
            self.disable_sliders()
        self.figure_3d.canvas.draw()
        plt.show()

    def on_pick_3d(self, event):
        self.ffd.on_pick(event)

    def disable_sliders(self):
        self.x_slide.configure(state="disabled")
        self.y_slide.configure(state="disabled")
        self.z_slide.configure(state="disabled")

    def enable_sliders(self):
        self.x_slide.configure(state="normal")
        self.y_slide.configure(state="normal")
        self.z_slide.configure(state="normal")

    def reset_sliders(self):
        self.x_slide.set(0)
        self.y_slide.set(0)
        self.z_slide.set(0)


if __name__ == "__main__":
    root = tk.Tk()
    app = MainApplication(root)
    app.start()
