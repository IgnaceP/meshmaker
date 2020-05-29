"""
Function to fill gaps(values '-inf') in a 2D Numpy Array

author: Ignace Pelckmans
                (University of Antwerp, Belgium)
"""

import numpy as np
from scipy.interpolate import griddata

def interpNumpy(arr, method = 'linear'):
    """
    Function fill gaps (values -inf) in a 2D Numpy array

    example:

        2    2    2    2    3    3    3         2    2    2    2    3    3    3
        2    2    2    2    3    3    3         2    2    2    2    3    3    3
      -inf -inf -inf -inf -inf -inf -inf  ==>   1    1    1    1   1.5  1.5  1.5
        0    0    0    0    0    0    0         0    0    0    0    0    0    0
        0    0    0    0    0    0    0         0    0    0    0    0    0    0

    Args:
        arr: (Required) Numpy array with dimensions m x n with gaps indicated with float('-inf')
        method: (Optional, defaults to 'linear') String with interpolation method, options: ‘linear’, ‘nearest’, ‘cubic’

    Returns:
        Numpy array with dimensions m x n with interpolated gaps
    """

    # get array dimensions
    rows, cols = np.shape(arr)

    # get array values as a flat array
    Z = arr.flatten()

    # meshgrid over the entire rasterized polygon array
    x = np.arange(arr.shape[1])
    y = np.arange(arr.shape[0])
    X, Y = np.meshgrid(x,y)

    X, Y = X.flatten(), Y.flatten()

    XY = np.zeros([len(Y), 2])
    XY[:,0] = X; XY[:,1] = Y

    # interpolate over the entire array
    Ti = griddata((X[Z > float('-inf')], Y[Z > float('-inf')]), Z[Z > float('-inf')], (X[Z == float('-inf')], Y[Z == float('-inf')]),
        method='linear')

    # replace -inf in flat array values
    Z[Z == float('-inf')] = Ti

    # Reshape the flatted array
    arr_interp = np.reshape(Z, [rows, cols])

    return arr_interp
