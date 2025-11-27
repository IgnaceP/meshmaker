"""
Script to declare function to create a background field, based on a geotiff file

author: Ignace Pelckmans
                    (University of Antwerp, Belgium)
"""

import numpy as np
from pathos.multiprocessing import ProcessingPool as PathosPool
import os
import shutil
import matplotlib.pyplot as plt
import matplotlib
from np2bgfield import np2BackgroundField

from osgeo import gdal, osr

def geotiff2Np(fn, return_projection = False, return_TL_coors = False, return_resolution = False):
    """
    Function to read a geotiff file and convert it to a Numpy array

    Args:
        fn: (Required) directory path string to indicate the filename of the raster
        return_projection: (Optional, defaults to False) True/False to indicate returning the epsg code
        return_TL_coors: (Optional, defaults to False) True/False to indicate returning the coordinates of the Top left corner
        return_resolution: (Optional, defaults to False) True/False to indicate returning the reoslution
    """

    # read file
    ds = gdal.Open(fn)
    # retreive epsg code
    proj = osr.SpatialReference(wkt=ds.GetProjection())
    epsg = proj.GetAttrValue('AUTHORITY',1)
    # retreive the coordinates of TL corner and resolution
    TL_x, res, _, TL_y, _, _ = ds.GetGeoTransform()
    # assume there is only one band
    band = ds.GetRasterBand(1)
    # read band as Numpy array
    arr = band.ReadAsArray()

    if return_projection == False and return_TL_coors == False and return_resolution == False:
        return arr
    else:
        ret = [arr]
        if return_projection: ret.append(epsg)
        if return_TL_coors: ret.append(TL_x);ret.append(TL_y)
        if return_resolution: ret.append(res)

        return ret

def geotiff2BackgroundField(fn, outputfile, multiplier = 1,  minres = 25, maxres = float('inf'), zerovalues = 250, flip_array = False):
    """
    Function to write a field file to accompagny a .geo file to determine resolution for GMSH based on a geotif

    author: Ignace Pelckmans
                        (University of Antwerp, Belgium)

    Args:
        fn: (Required) directory path string to indicate the filename of the raster
        outputfile: (Required) directory to store the newly generated file
        multiplier: (Optional, defaults 1) value to multiply the original array with
        minres: (Optional, defaults to 25) minimum mesh size
        maxres: (Optional, defaults to inf) maximum mesh size
        zerovalues: (Optional, defaults to 250) the mesh size in areas with value 0
    """
    xy, TL_x, TL_y, res = geotiff2Np(fn, return_TL_coors = True, return_resolution = True)
    rows, cols = xy.shape
    BL_x = TL_x
    BL_y = TL_y - res*rows

    if flip_array:
        xy = np.flip(xy, axis = 0)

    np2BackgroundField(xy, outputfile, multiplier = multiplier, minres = minres, maxres = maxres, zerovalues = zerovalues, zeropoint = [TL_x,TL_y], xres = res, yres = res)
    print('Succesfully created a backgroundfield file in %s' % outputfile)

#geotiff2BackgroundField('/media/ignace/LaCie/TELEMAC/Churute_setup/gmsh/background_field.tif', '/media/ignace/LaCie/TELEMAC/Churute_setup/gmsh/background_field.bg')
