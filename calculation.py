import numpy as np


def vectorProjection(vector, base_vector):
    """returns the projection value of vector on base_vector"""
    return np.dot(vector, base_vector) / np.linalg.norm(base_vector)

