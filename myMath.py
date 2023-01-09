# Custom math functions for compression algorithms #
import numpy as np
import math

def invertMatrix(matrix):
    """
    Takes in matrix (multidimentional python list)

    Outputs matrix (multidimentional python list)
    """
    mat = [[0]*len(matrix) for _ in range(len(matrix[0]))] #initialize
    for row in range(len(matrix)):
        for col in range(len(matrix[row])):
            mat[col][row] = matrix[row][col]
    
    return mat