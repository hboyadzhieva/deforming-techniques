import matplotlib.pyplot as plt

from deformator import GridDeformator2D, FreeFormDeformator
from model import Model2D, Model3D
from point import Point2D, Point3D

ax = None


def main():
    test_model = Model2D(Point2D(2.5, 4))
    print(test_model)
    test_deformator = GridDeformator2D(test_model, 5, 5)
    print(test_deformator)

    test_model_3d = Model3D(Point3D(1.5, 0, 2.5), Point3D(0, 1.5, 1.5))
    test_deformator_3d = FreeFormDeformator(test_model_3d)
    print(test_deformator_3d)
    # fig = plt.figure(figsize=(4, 4))
    # global ax
    # ax = fig.add_subplot(111, projection='3d')
    # ax.scatter(2, 3, 4)
    # fig.canvas.callbacks.connect('button_press_event', onmousedown)
    # plt.show()


def onmousedown(event):
    global ax
    result = ax.format_coord(event.xdata, event.ydata);
    result.split(',')
    print(result)


# Press the green button in the gutter to run the script.
if __name__ == '__main__':
    main()
