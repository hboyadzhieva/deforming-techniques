from deformator import GridDeformator2D, FreeFormDeformator
from model import Model2D, Model3D
from point import Point2D, Point3D

ax = None


def main():
    test_model = Model2D(Point2D(2.5, 4))
    test_deformator = GridDeformator2D(test_model, 5, 5)
    test_deformator.start()

    test_model_3d = Model3D(Point3D(1, 0, 2), Point3D(3, 2, 5))
    test_deformator_3d = FreeFormDeformator(test_model_3d)
    test_deformator_3d.start()


def onmousedown(event):
    global ax
    result = ax.format_coord(event.xdata, event.ydata)
    result.split(',')
    print(result)


# Press the green button in the gutter to run the script.
if __name__ == '__main__':
    main()
