import tkinter as tk

import matplotlib.pyplot as plt
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg

from model import Model2D, Model3D
from deformator import GridDeformator2D, FreeFormDeformator


def display_option():
    global canvas, model_2d_option_1, model_2d_option_2, \
        model_3d_option_1, model_3d_option_2, x_slide, y_slide, z_slide, show_grid
    canvas.get_tk_widget().pack_forget()
    if selected_def_option.get() == 1:
        model_3d_option_1.pack_forget()
        model_3d_option_2.pack_forget()
        x_slide.pack_forget()
        y_slide.pack_forget()
        z_slide.pack_forget()
        model_2d_option_1.pack()
        model_2d_option_2.pack()
        show_grid.pack()
    elif selected_def_option.get() == 2:
        model_2d_option_1.pack_forget()
        model_2d_option_2.pack_forget()
        model_3d_option_1.pack()
        model_3d_option_2.pack()
        x_slide.pack()
        y_slide.pack()
        z_slide.pack()
        show_grid.pack()


def start_2d_deformator_circle():
    global canvas, ax_2d, def_2d
    canvas.get_tk_widget().pack_forget()
    ax_2d.cla()
    canvas = FigureCanvasTkAgg(figure_2d, root)
    canvas.get_tk_widget().pack()
    circle = Model2D.createCircle(2)
    def_2d = GridDeformator2D(circle)
    def_2d.plot_model(ax_2d)
    def_2d.plot_grid(ax_2d)
    plt.show()


def start_2d_deformator_triangle():
    global canvas, def_2d, ax_2d
    canvas.get_tk_widget().pack_forget()
    ax_2d.cla()
    canvas = FigureCanvasTkAgg(figure_2d, root)
    canvas.get_tk_widget().pack()
    sin_function = Model2D.createSin()
    def_2d = GridDeformator2D(sin_function)
    def_2d.plot_model(ax_2d)
    def_2d.plot_grid(ax_2d)
    plt.show()


def start_3d_deformator_sphere():
    global canvas, def_3d, ax_3d
    canvas.get_tk_widget().pack_forget()
    canvas = FigureCanvasTkAgg(figure_3d, root)
    canvas.get_tk_widget().pack()
    ax_3d.cla()
    sphere = Model3D.sphere()
    def_3d = FreeFormDeformator(sphere)
    def_3d.plot_model(ax_3d)
    def_3d.plot_grid(ax_3d)
    plt.show()


def start_3d_deformator_torus():
    global canvas, def_3d, ax_3d
    canvas.get_tk_widget().pack_forget()
    canvas = FigureCanvasTkAgg(figure_3d, root)
    canvas.get_tk_widget().pack()
    ax_3d.cla()
    torus = Model3D.torus()
    def_3d = FreeFormDeformator(torus)
    def_3d.plot_model(ax_3d)
    def_3d.plot_grid(ax_3d)
    plt.show()


def on_mouse_down_2d(event):
    global def_2d
    def_2d.on_mouse_down(event)


def on_mouse_up_2d(event):
    global def_2d
    def_2d.on_mouse_up(event)


def on_mouse_motion_2d(event):
    global ax_2d, def_2d, figure_2d
    def_2d.on_mouse_motion(event)
    ax_2d.cla()
    def_2d.plot_model(ax_2d)
    def_2d.plot_grid(ax_2d)
    figure_2d.canvas.draw()
    plt.show()


def on_mouse_down_3d(event):
    global def_3d, ax_3d, x_slide, show_grid_selected
    def_3d.on_mouse_down(event, ax_3d)
    reset_sliders()
    if def_3d._vertex_selected is not None:
        ax_3d.cla()
        def_3d.add_arrows(ax_3d)
        def_3d.plot_model(ax_3d)
        if show_grid_selected.get() == 1:
            def_3d.plot_grid(ax_3d)
        x_slide.configure(state="normal")
        y_slide.configure(state="normal")
        z_slide.configure(state="normal")
    else:
        ax_3d.cla()
        def_3d.plot_model(ax_3d)
        if show_grid_selected.get() == 1:
            def_3d.plot_grid(ax_3d)
        x_slide.configure(state="disabled")
        y_slide.configure(state="disabled")
        z_slide.configure(state="disabled")
    figure_3d.canvas.draw()
    plt.show()


def on_pick_3d(event):
    global def_3d
    def_3d.on_pick(event)


def change_x_pos(event):
    change_vertex_position(1)


def change_y_pos(event):
    change_vertex_position(2)


def change_z_pos(event):
    change_vertex_position(3)


def change_vertex_position(coord_axis):
    ax_3d.cla()
    if def_3d._vertex_selected is not None and def_3d._original_vertex_selected is not None:
        if coord_axis == 1:
            def_3d._vertex_selected.x = def_3d._original_vertex_selected.x + dx.get()
        if coord_axis == 2:
            def_3d._vertex_selected.y = def_3d._original_vertex_selected.y + dy.get()
        if coord_axis == 3:
            def_3d._vertex_selected.z = def_3d._original_vertex_selected.z + dz.get()
        def_3d.add_arrows(ax_3d)
    def_3d.plot_model(ax_3d)
    if show_grid_selected.get() == 1:
        def_3d.plot_grid(ax_3d)
    figure_3d.canvas.draw()
    plt.show()


def toggle_grid():
    global show_grid_selected, selected_def_option, def_3d, def_3d, ax_3d, ax_2d
    if selected_def_option.get() == 1:
        ax_2d.cla()
        def_2d.plot_model(ax_2d)
        if show_grid_selected.get() == 1:
            def_2d.plot_grid(ax_2d)
        figure_2d.canvas.draw()
        plt.show()
    elif selected_def_option.get() == 2:
        ax_3d.cla()
        def_3d.plot_model(ax_3d)
        if show_grid_selected.get() == 1:
            def_3d.plot_grid(ax_3d)
        figure_3d.canvas.draw()
        plt.show()


def reset_sliders():
    global x_slide, y_slide, z_slide
    x_slide.set(0)
    y_slide.set(0)
    z_slide.set(0)


root = tk.Tk()
root.geometry('%dx%d+%d+%d' % (1000, 800, 0, 0))
root.title('Object Deformation Algorithms Simulation')

selected_def_option = tk.IntVar()
def_options = (("2D Grid Deformator", 1), ("Free Form Deformator", 2))
for option in def_options:
    r = tk.Radiobutton(root, text=option[0], value=option[1], variable=selected_def_option, command=display_option)
    r.pack(padx=5, pady=5)

figure_2d = plt.Figure(figsize=(6, 5), dpi=100)
ax_2d = figure_2d.add_subplot(111)
def_2d = None
selected_model_2d_option = tk.IntVar()
model_2d_option_1 = tk.Radiobutton(root, text="Circle", value=1, variable=selected_model_2d_option,
                                   command=start_2d_deformator_circle)
model_2d_option_2 = tk.Radiobutton(root, text="Triangle", value=2, variable=selected_model_2d_option,
                                   command=start_2d_deformator_triangle)
figure_2d.canvas.mpl_connect('button_press_event', on_mouse_down_2d)
figure_2d.canvas.mpl_connect('button_release_event', on_mouse_up_2d)
figure_2d.canvas.mpl_connect('motion_notify_event', on_mouse_motion_2d)

figure_3d = plt.Figure(figsize=(6, 5), dpi=100)
ax_3d = figure_3d.add_subplot(111, projection='3d')
def_3d = None
selected_model_3d_option = tk.IntVar()
model_3d_option_1 = tk.Radiobutton(root, text="Sphere", value=1, variable=selected_model_3d_option,
                                   command=start_3d_deformator_sphere)
model_3d_option_2 = tk.Radiobutton(root, text="Torus", value=2, variable=selected_model_3d_option,
                                   command=start_3d_deformator_torus)

dx, dy, dz = tk.DoubleVar(), tk.DoubleVar(), tk.DoubleVar()
x_slide = tk.Scale(root, from_=-2, to=2, resolution=0.01, orient=tk.HORIZONTAL, variable=dx,
                   command=change_x_pos, fg="red", label="x")
y_slide = tk.Scale(root, from_=-2, to=2, resolution=0.01, orient=tk.HORIZONTAL, variable=dy,
                   command=change_y_pos, fg="blue", label="y")
z_slide = tk.Scale(root, from_=-2, to=2, resolution=0.01, orient=tk.HORIZONTAL, variable=dz,
                   command=change_z_pos, fg="green", label="z")
x_slide.configure(state="disabled")
y_slide.configure(state="disabled")
z_slide.configure(state="disabled")

figure_3d.canvas.mpl_connect('button_press_event', on_mouse_down_3d)
figure_3d.canvas.mpl_connect('pick_event', on_pick_3d)

show_grid_selected = tk.IntVar(value=1)
show_grid = tk.Checkbutton(root, text='Show gridlines', variable=show_grid_selected, onvalue=1, offvalue=0,
                           command=toggle_grid)

canvas = FigureCanvasTkAgg(figure_2d, root)
root.mainloop()
