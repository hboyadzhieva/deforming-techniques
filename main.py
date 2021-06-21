import matplotlib.pyplot as plt

from deformator import GridDeformator2D
from model import Model2D
from point import Point2D

ax = None


def main():
    test_model = Model2D(Point2D(2.5, 4))
    print(test_model)
    test_deformator = GridDeformator2D(test_model, 5, 5)
    print(test_deformator)
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
