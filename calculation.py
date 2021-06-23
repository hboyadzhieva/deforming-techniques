import math

import numpy as np


def vector_projection(vector, base_vector):
    """Calculate the projection of vector on base_vector

    Parameters:
    ----------
    vector: numpy.array
        vector to project
    base_vector: numpy.array\
        base vector projection is made on

    Returns:
    -------
    float:
        numeric value representing length of projection of vector on base_vector
    """
    return np.dot(vector, base_vector) / np.linalg.norm(base_vector)


def trilinear_interpolant(A, B, position_vector, base_vector):
    """Calculate (AxB).(position_vector)/(AxB).(base_vector)

    Parameters
    ----------
    A: numpy.array
    B: numpy.array
        AxB creates orthogonal vector to A and B
    base_vector: numpy.array
        the vector we are projecting onto
    position_vector: numpy.array
        the vector being projected

    Returns
    -------
    float
        numeric value representing the projection of position_vector on base_vector
    """
    return np.dot(np.cross(A, B), position_vector) / np.dot(np.cross(A, B), base_vector)


def bilinear_interpolation(P00, P10, P11, P01, h, v):
    """Calculate the position of a point on a curve defined by bilinear interploation using four countrol points

    Parameters
    -----------
    P00: numpy.array
        vector representing left vertex on lower line
    P10: numpy.array
        vector representing right vertex on lower line
    P01: numpy.array
        vector representing left vertex on upper line
    P11: numpy.array
        vector representing right vertex on upper line
    h: float
        horizontal interpolant
    v: float
        vertical interpolant

    Returns
    --------
        array representing the 2D result vector
    """
    U = (1 - h) * P00 + h * P10
    V = (1 - h) * P01 + h * P11
    A = (1 - v) * U + v * V
    return A


def combination(n, i):
    """Number of possible arrangements by selecting only a few elements from a set

    Parameters
    ---------
    n: int
        all elements in set
    i: int
        the number of selected objects

    Returns
    --------
    int:
        the number of all possible combinations
    """
    return math.factorial(n) / (math.factorial(n - i) * math.factorial(i))


def bernstain(n, i, t):
    """ Value of i-th Bernstain polynomial of degree n at t"""
    return combination(n, i) * ((1 - t) ** (n - i)) * (t ** i)


def bezier_volume(control_grid, s, t, u):
    result = np.zeros(3)
    for i in range(0, control_grid.count_x_points):
        for j in range(0, control_grid.count_y_points):
            for k in range(0, control_grid.count_z_points):
                result += control_grid.control_points[i][j][k].to_numpy_array() * \
                          bernstain(control_grid.count_x_points - 1, i, s) * \
                          bernstain(control_grid.count_y_points - 1, j, t) * \
                          bernstain(control_grid.count_z_points - 1, k, u)
    return result
