from deformator import GridDeformator2D, FreeFormDeformator
from model import Model2D, Model3D
from point import Point2D, Point3D

ax = None


def main():
    test_model = Model2D.createCircle(2)
    test_deformator = GridDeformator2D(test_model)
    # test_deformator.start()

    test_model_3d = Model3D.torus()
    test_deformator_3d = FreeFormDeformator(test_model_3d,offset_grid_model=0)
    test_deformator_3d.start()



if __name__ == '__main__':
    main()
