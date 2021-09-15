"""Util module to read files from the file system"""
import numpy as np


def read_obj(filename):
    """
    Reads .obj file from file system and returns vertices and faces data.

    :param filename: absolute path to .obj file
    :return: array of vertices of the object (numpy.array),
             array of faces (triangles) of the object (numpy.array)

    *********************************************************************************************************************
    Title: Beyond data scientist: 3d plots in Python with examples
    Author: Yuchen Zhong
    Date: Jul 18, 2019
    Availability: https://yzhong-cs.medium.com/beyond-data-scientist-3d-plots-in-python-with-examples-2a8bd7aa654b
    *********************************************************************************************************************
    """
    triangles = []
    vertices = []
    with open(filename) as file:
        for line in file:
            components = line.strip(' \n').split(' ')
            if components[0] == "f":  # face data
                # e.g. "f 1/1/1/ 2/2/2 3/3/3 4/4/4 ..."
                indices = list(map(lambda c: int(c.split('/')[0]) - 1, components[1:]))
                for i in range(0, len(indices) - 2):
                    triangles.append(indices[i: i + 3])
            elif components[0] == "v":  # vertex data
                # e.g. "v  30.2180 89.5757 -76.8089"
                vertex = list(map(lambda c: float(c), components[1:]))
                vertices.append(vertex)
    return np.array(vertices), np.array(triangles)
