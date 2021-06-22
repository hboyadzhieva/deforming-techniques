import numpy as np


def vector_projection(vector, base_vector):
    """
    Returns the projection value of vector on base_vector
    Parameters:
        vector: numpy array, vector to project on other
        base_vector: numpy array, base vector projection is made on
    """
    return np.dot(vector, base_vector) / np.linalg.norm(base_vector)


def trilinear_interpolant(A, B, position_vector, base_vector):
    """
    Return (AxB).(position_vector)/(AxB).(base_vector)
    Parameters:
        A: numpy array, AxB creates orthogonal vector to A and B
        B: numpy array, AxB creates orthogonal vector to A and B
        base_vector: numpy array, the vector we are projecting onto
        position_vector: numpy array, the vector being projected
    """
    return np.dot(np.cross(A, B), position_vector) / np.dot(np.cross(A, B), base_vector)
