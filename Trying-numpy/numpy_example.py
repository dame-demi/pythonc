# numpy_example.py

import numpy as np
def numpy_function(matrix):
    """ Calculate the sum of all elements in the given matrix """
    np_matrix = np.array(matrix)
    result = np.sum(np_matrix)
    print(f"The sum of all elements in the matrix is: {result}")
    return result


