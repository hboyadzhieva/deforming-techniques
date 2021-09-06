"""Helper module for calculations in 2D and 3D"""
import math
import numpy as np


def distance_2d(vector1, vector2):
    """
    Calculate the distance between 2 points in 2D
    @:param vector1 (point.Point2D), vector2 (point.Point2D)
    @:returns distance between vector1 and vector2 (float)
    """
    return math.sqrt((vector2.x - vector1.x) ** 2 + (vector2.y - vector1.y) ** 2)


def distance_3d(vector1, vector2):
    """
    Calculate the distance between 2 points in 3D
    @:param vector1 (point.Point3D), vector2 (point.Point3D)
    @:returns distance between vector1 and vector2 (float)
    """
    return math.sqrt((vector2.x - vector1.x) ** 2 + (vector2.y - vector1.y) ** 2 + (vector2.z - vector1.z) ** 2)


def scalar_projection(vector, base_vector):
    """
    Calculate the projection of vector on base_vector
    @:param vector (numpy.array), base_vector (numpy.array)
    @:returns numeric value representing length of projection of vector on base_vector (float)
    """
    return np.dot(vector, base_vector) / np.linalg.norm(base_vector)


def trilinear_interpolant(A, B, position_vector, C):
    """
    Calculate the local coordinate of position_vector in accordance to C,
    where A, B and C are bases of new coordinate system
    @:param A (numpy.array), B (numpy.array), C (numpy.array), position_vector (numpy.array)
    @:returns numeric value representing the projection of position_vector on C (float)
    """
    return np.dot(np.cross(A, B), position_vector) / np.dot(np.cross(A, B), C)


def bilinear_interpolation(P00, P10, P11, P01, h, v):
    """
    Calculate the position of a point on a curve defined by bilinear interpolation using four countrol points
    @:param P00: lower left vertex (numpy.array), P10: lower right vertex (numpy.array)
    @:param P01: upper left vertex (numpy.array), P11: upper right vertex (numpy.array)
    @:param h: horizontal interpolant (float), v: vertical interpolant (float)
    @:returns the 2D vector representing the new position of the point on the curve (numpy.array)
    """
    U = (1 - h) * P00 + h * P10
    V = (1 - h) * P01 + h * P11
    A = (1 - v) * U + v * V
    return A


def combination(n, i):
    """Number of possible arrangements by selecting i elements from a total of n"""
    return math.factorial(n) / (math.factorial(n - i) * math.factorial(i))


def bernstain(n, i, t):
    """Value of i-th Bernstain polynomial of degree n at t"""
    return combination(n, i) * ((1 - t) ** (n - i)) * (t ** i)


def point_position_in_bezier_volume(control_grid, s, t, u):
    """
    Calculates position of a point given its trilinear interpolants (s,t,u)
    and the control grid that is defining the bezier volume
    """
    result = np.zeros(3)
    for i in range(0, control_grid.count_x_points):
        for j in range(0, control_grid.count_y_points):
            for k in range(0, control_grid.count_z_points):
                result += control_grid.control_points[i][j][k].to_numpy_array() * \
                          bernstain(control_grid.count_x_points - 1, i, s) * \
                          bernstain(control_grid.count_y_points - 1, j, t) * \
                          bernstain(control_grid.count_z_points - 1, k, u)
    return result
