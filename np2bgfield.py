"""
Script to declare function to create a background field, based on a Numpy 2D array

author: Ignace Pelckmans
                    (University of Antwerp, Belgium)
"""

import numpy as np
from pathos.multiprocessing import ProcessingPool as PathosPool
import os
import shutil
import matplotlib.pyplot as plt
import matplotlib
from scipy.ndimage import morphology

from interpolate2DNumpy import *


def np2BackgroundField(xy, outputfile, multiplier = 1,  minres = 25, maxres = float('inf'), zerovalues = 250, zeropoint = [0,0],
    xres = 1, yres = 1, plot_arr = False, buffer_interpolation = False):
    """
    Function to write a field file to accompagny a .geo file to determine resolution for GMSH

    author: Ignace Pelckmans
                        (University of Antwerp, Belgium)

    Args:
        xy: (Required) Numpy Array with dimension m x n with the values of mesh size
        outputfile: (Required) directory to store the newly generated file
        multiplier: (Optional, defaults 1) value to multiply the original array with
        minres: (Optional, defaults to 25) minimum mesh size
        maxres: (Optional, defaults to inf) maximum mesh size
        zerovalues: (Optional, defaults to 250) the mesh size in areas with value 0
        zeropoint: (Optional, defaults to (0,0)) coordinate pair of the left bottom corner of the numpy image
        xres: (Optional, defaults to 1) cell width
        yres: (Optional, defaults to 1) cell height
        plot_arr: (Optional, defaults to False) False/file directory to store a plot of the original raster
        buffer_interpolation: (Optional, defaults to False) False/float to indicate the size of buffer around the channels to interpolate in meters
    """

    print('Translating the numpy array into a background file...')
    # temporary directory
    if os.path.isdir('./tmp/'): shutil.rmtree('./tmp/')
    os.mkdir('./tmp/')

    # copy original
    xy_orig = xy.copy()

    # apply multiplier:
    xy = xy*multiplier

    # assign zerovalues
    xy[xy_orig == 0] = zerovalues
    # if buffer is indicated:
    if buffer_interpolation:
        channels = xy_orig > 0
        dist = morphology.distance_transform_edt(channels==0)
        dist[dist == 0] = float('inf')
        xy[dist <= buffer_interpolation/xres] = float('-inf')
        xy = interpNumpy(xy, method = 'cubic')

    # replace nans with the zerovalue
    xy[np.isnan(xy)] = zerovalues
    # replace all values smaller than the minimum cell size with the minimum resolution
    xy[xy < minres] = minres
    # replace all values larger than the maximum cell size with the maximum resolution
    xy[xy > maxres] = maxres

    # if indicated plot a figure with the field
    if plot_arr:
        f = plt.figure(figsize = [12,12])
        xy_adapted = xy.copy()
        xy_adapted[xy_orig != 0] = 0
        ims = plt.imshow(xy, norm=matplotlib.colors.LogNorm())
        plt.colorbar(ims)
        f.savefig(plot_arr)

    # open a file to store the header and allow writing
    f = open('./tmp/header_file.txt', "w+")
    # write the first line which are a coordinate triple with the location of the left bottom
    f.write('%f %f %f\n' % (zeropoint[0] + xres / 2, zeropoint[1] + yres / 2, 0))
    # write the second line which are the resolutions
    f.write('%f %f %f\n' % (xres, yres, 1))
    # write the thrid line which are the total number of cells in each the x-, y- and z-direction
    f.write('%d %d %d\n' % (np.shape(xy)[1], np.shape(xy)[0], 1))
    # close the header line
    f.close()

    # flatten the 2D Numpy array
    xy = np.flip(xy.T, axis=1).flatten()


    # create a function to write the values per line without a header to enable parallelizing
    def writeFile(n, lis):
        # open a file to store the values
        f = open('./tmp/file_%s.txt' % (str(n)), "w+")
        for t in range(len(lis)):
            # select the value
            v = lis[t]
            # it is a 2D field, so only one value per line
            f.write('%f\n' % (v))
        # close the open file
        f.close()
        # return the name of the newly written file
        return './tmp/file_%s.txt' % (str(n))

    # divide all values of the original array into 8 blocks
    n, n_rest = len(xy) // 8, len(xy) % 8
    N, N[0], N[-1] = [n for i in range(9)] * np.arange(9), 0, len(xy)
    xy_tiles = [xy[N[i]:N[i + 1]] for i in range(8)]

    print('writing the separate txt files...')
    # parallelize the function using PathosPool
    pool = PathosPool(8)
    P = pool.map(writeFile, np.arange(8), xy_tiles)
    # add the header file as the first txt file in the list
    P.insert(0,'./tmp/header_file.txt')

    print('stitching them together...')
    # stitch all txt files to each other, write in binary mode
    with open(outputfile, 'wb') as wfd:
        # open each txt file and copy its content into the final txt file
        for f in P:
            with open(f, 'rb') as fd:
                shutil.copyfileobj(fd, wfd)

    # remove the tmp directory
    shutil.rmtree('./tmp/')
    print('Succes!')
