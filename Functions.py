""" Script to store all functions"""

import numpy as np
from osgeo.gdalconst import *
from osgeo import ogr, gdal
import subprocess
import simplekml
import sys
import scipy
import multiprocessing as mp
from scipy import ndimage
import os
import cv2
import time
import multiprocessing as mp
from pathos.multiprocessing import ProcessingPool as PathosPool
import functools
from scipy.spatial import Voronoi, voronoi_plot_2d
import matplotlib.pyplot as plt
from shapely.geometry import Polygon as shPol
from shapely.geometry.polygon import LinearRing as shLinRing
from matplotlib.collections import PatchCollection
from matplotlib.patches import Polygon
from matplotlib.path import Path
from shapely.geometry import Polygon as shPol
from shapely.geometry import MultiPolygon as shMPol
from shapely.geometry import Point as shPoint
from osgeo import ogr
import json
from shapely.geometry import shape
import meshio
import numpy as np
import pygmsh
import re
import utm
import shutil
from osgeo import osr

def distanceMapCenter(buffer_size):
    """ function to create a raster with distances to the central cell,
        rounded, only input is the size of the buffer around the middle cell """
    
    rows = buffer_size*2+1
    dist = np.zeros([rows, rows])
    mid = buffer_size # buffer size is alway the index of the middle cell
    
    for i in range(0,rows):
        for j in range(0,rows):
            dist[i,j] = np.sqrt((i-mid)**2+(j-mid)**2)
            dist = np.round(dist)
        
    return(dist)

def rasterize(InputVector,RefImage,OutputImage):
    """ function to rasterize vectorial data, works with polygons, points and lines
    all corresponding cells get value 1 """

    #InputVector = 'Roads.shp'
    #OutputImage = 'Result.tif'
    
    #RefImage = 'DEM.tif'
    
    gdalformat = 'GTiff'
    datatype = gdal.GDT_Byte
    burnVal = 1 #value for the output image pixels
    ##########################################################
    # Get projection info from reference image
    Image = gdal.Open(RefImage, gdal.GA_ReadOnly)
    
    # Open Shapefile
    Shapefile = ogr.Open(InputVector)
    Shapefile_layer = Shapefile.GetLayer()
    
    # Rasterise
    print("Rasterising shapefile...")
    Output = gdal.GetDriverByName(gdalformat).Create(OutputImage, Image.RasterXSize, Image.RasterYSize, 1, datatype, options=['COMPRESS=DEFLATE'])
    Output.SetProjection(Image.GetProjectionRef())
    Output.SetGeoTransform(Image.GetGeoTransform()) 
    
    # Write data to band 1
    Band = Output.GetRasterBand(1)
    Band.SetNoDataValue(0)
    gdal.RasterizeLayer(Output, [1], Shapefile_layer, burn_values=[burnVal])
    
    # Close datasets
    Band = None
    Output = None
    Image = None
    Shapefile = None
    
    # Build image overviews
    subprocess.call("gdaladdo --config COMPRESS_OVERVIEW DEFLATE "+OutputImage+" 2 4 8 16 32 64", shell=True)
    print("Done.")
    
    return(OutputImage)

def localEdgeRasterize(XY, generalize_factor=5, get_indices = False):
    """ Function to locally rasterize a set of coordinates
    @params:
    XY                  -       Required: List of coordinate pairs
    generalize_factor   -       Optional: Factor to determine the degree of generalization"""

    XY = addVerticesToPolygon(XY, generalize_factor, shapely= False)
    XY_upd = []
    for x, y in XY:
        x = int(x / generalize_factor)
        y = int(y / generalize_factor)
        XY_upd.append([x, y])

    XY_np = np.asarray(XY_upd)
    x_min = np.min(XY_np[:, 0])
    x_max = np.max(XY_np[:, 0])
    y_min = np.min(XY_np[:, 1])
    y_max = np.max(XY_np[:, 1])
    cols = int(x_max - x_min + 2)
    rows = int(y_max - y_min + 2)

    ras = np.zeros([rows, cols])

    XY_in = []
    for c, r in XY_upd:
        c = int(c - x_min + 1)
        r = int(y_max - r + 1)
        XY_in.append([r, c])
        ras[r, c] = 1

    if get_indices:
        return ras, XY_in
    else:
        return ras

def localRasterize(XY_exter_L, XY_inner_L = [], generalize_factor = 1, get_indices = False, xmin = [], xmax = [], ymin = [], ymax = [], printing = False ):
    """ Function to locally rasterize a set of coordinates
    @params:
    XY_exter_L          -       Required: List of coordinate pairs of outer bounds
    XYs_inner           -       Optional: List of coordinate pairs lists of inner bounds

   """

    t_start = time.time()

    if XY_exter_L == []:
        if type(xmin) != list and type(xmax) != list and type(ymin) != list and type(ymax) != list:
            print('Woeps Empty Polygon, zero array is returned.')
            return np.zeros([ymax-ymin, xmax-xmin])
        else:
            print('Woeps Empty Polygon, 0 is returned.')
            return(0)
    else:
        # make sure it can accept both a single coordiante pair list or a list of lists
        if len(XY_exter_L[0]) == 2:
            XY_exter_L = [XY_exter_L]

        if len(XY_inner_L) > 0:
            if len(XY_inner_L[0]) == 2:
                XY_inner_L = [XY_inner_L]

        XY_exter_L_upd = []
        XY_inner_L_upd = []

        x_mins = []
        y_mins = []
        x_maxs = []
        y_maxs = []

        for XY_exter in XY_exter_L:
            XY_exter_upd = []
            for x, y in XY_exter:
                x = int(x / generalize_factor)
                y = int(y / generalize_factor)
                XY_exter_upd.append([x, y])
            XY_exter_L_upd.append(XY_exter_upd)

            XY_exter_upd_np = np.asarray(XY_exter_upd)
            x_mins.append(np.min(XY_exter_upd_np[:, 0]))
            x_maxs.append(np.max(XY_exter_upd_np[:, 0]))
            y_mins.append(np.min(XY_exter_upd_np[:, 1]))
            y_maxs.append(np.max(XY_exter_upd_np[:, 1]))

        if len(XY_inner_L) > 0:
            for XY_inner in XY_inner_L:
                XYs_inner_upd = []
                for i in XY_inner:
                    XYs_inner_upd.append([])
                    for xy in i:
                        x = xy[0]
                        y = xy[1]
                        x = int(x / generalize_factor)
                        y = int(y / generalize_factor)
                        XYs_inner_upd[-1].append([x, y])
                XY_inner_L_upd.append(XYs_inner_upd)

        if type(xmax) == list: x_max = np.max(x_maxs)
        else: x_max = xmax
        if type(ymax) == list: y_max = np.max(y_maxs)
        else: y_max = ymax
        if type(xmin) == list: x_min = np.min(x_mins)
        else: x_min = xmin
        if type(ymin) == list: y_min = np.min(y_mins)
        else: y_min = ymin
        cols = int(x_max - x_min)
        rows = int(y_max - y_min)

        #print('xmax: ', x_max, 'and xmin: ', x_min)
        #print('ymax: ', y_max, 'and ymin: ', y_min)


        # Create vertex coordinates for each grid cell...
        # (<0,0> is at the top left of the grid in this system)
        print('preparing the mesh grid...')
        x_ind, y_ind = np.meshgrid(np.arange(x_min, x_max), np.arange(y_min, y_max))
        x, y = np.meshgrid(np.arange(cols), np.arange(rows))
        x = x.flatten()+x_min
        y = y.flatten()+y_min
        points = np.vstack((x, y))
        points = points.T

        # outer bounds
        t = 0
        print('rasterizing...')
        for XY in XY_exter_L_upd:
            t += 1
            printProgressBar(t, len(XY_inner_L_upd)+1)
            path = Path(XY)
            grid = path.contains_points(points)
            grid = grid.reshape((rows, cols))
            if t == 1:
                grid_outer = grid
            grid_outer = grid|grid_outer


        # inner:
        t = 0
        if len(XY_inner_L) > 0:
            for XY in XY_inner_L_upd:
                    for I in XY:
                        if len(I) > 0:
                            t += 1
                            path = Path(I)
                            grid = path.contains_points(points)
                            grid = grid.reshape((rows, cols))
                            grid = np.invert(grid)
                            grid_outer *= grid

        grid_outer = np.flip(np.flip(grid_outer), axis = 1)

        if get_indices:
            return grid_outer.astype(int), x_ind[0,:], y_ind[:,0]
        else:
            return grid_outer.astype(int)

        tit = time.time() - t_start
        min = int(tit // 60)
        sec = int(tit % 60)
        print('It took %d minutes and %d seconds to rasterize this tile!' % (min, sec))

def localRasterizeInParallel(XY_exter_L, XY_inner_L = [], generalize_factor=1, get_indices = False):
    """ Function to locally rasterize a set of coordinates
    @params:
    XY_exter_L          -       Required: List of coordinate pairs of outer bounds
    XYs_inner           -       Optional: List of coordinate pairs lists of inner bounds
    generalize_factor   -       Optional: Factor to determine the degree of generalization"""

    # track time
    t_start = time.time()

    # make sure it can accept both a single coordiante pair list or a list of lists
    if len(XY_exter_L[0]) == 2:
        XY_exter_L = [XY_exter_L]

    if len(XY_inner_L) > 0:
        if len(XY_inner_L[0][0]) == 2:
            XY_inner_L = [XY_inner_L]
            inner_flag = True
    else:
        for i in range(len(XY_exter_L)):
            XY_inner_L.append([])
            inner_flag = False

    # generalize
    XY_exter_L_upd = []
    for pol in XY_exter_L:
        XY_exter_L_upd.append([])
        for x,y in pol:
            x = x/generalize_factor
            y= y/generalize_factor
            XY_exter_L_upd[-1].append((x,y))


    XY_exter_L = XY_exter_L_upd

    XY_inner_L_upd = []
    if len(XY_inner_L) > 0:
        for XY_inner in XY_inner_L:
            XYs_inner_upd = []
            for i in XY_inner:
                XYs_inner_upd.append([])
                for x, y in i:
                    x = x / generalize_factor
                    y = y / generalize_factor
                    XYs_inner_upd[-1].append((x, y))

            XY_inner_L_upd.append(XYs_inner_upd)

    XY_inner_L = XY_inner_L_upd

    # calculate the real mins and maxs
    x_mins, y_mins, x_maxs, y_maxs = [],[],[],[]
    for XY_exter in XY_exter_L:
        XY_exter_np = np.asarray(XY_exter)
        x_mins.append(np.min(XY_exter_np[:, 0]))
        x_maxs.append(np.max(XY_exter_np[:, 0]))
        y_mins.append(np.min(XY_exter_np[:, 1]))
        y_maxs.append(np.max(XY_exter_np[:, 1]))

    x_max, y_max, x_min, y_min = np.max(x_maxs), \
                                 np.max(y_maxs), \
                                 np.min(x_mins), \
                                 np.min(y_mins)

    # calculate dimensions
    cols, rows = int(x_max - x_min), \
                 int(y_max - y_min)

    print('xmax: ', x_max, 'and xmin: ', x_min) ; print('ymax: ', y_max, 'and ymin: ', y_min)

    # original indices
    xmax_orig = x_max*generalize_factor
    ymax_orig = y_max*generalize_factor
    xmin_orig = x_min*generalize_factor
    ymin_orig = y_min*generalize_factor

    x_coordinates = np.hstack((np.arange(xmin_orig, xmax_orig, (xmax_orig-xmin_orig)/cols),xmax_orig))
    y_coordinates = np.hstack((np.arange(ymin_orig, ymax_orig, (ymax_orig-ymin_orig)/rows),ymax_orig))

    # create list of shapely polygons
    Pols = [] ; t = 0
    for p in XY_exter_L:
        if len(XY_inner_L) > 0:
            inter = XY_inner_L[t]

        inter = [i for i in inter if len(i) >= 3]
        if len(inter) > 0: pol = shPol(p, inter)
        else: pol = shPol(p)
        Pols.append(pol)
        t += 1
    MultiPol = shMPol(Pols)

    # calculate tile dimensions
    width = cols // 4
    height = rows // 2
    print('width: ', cols , ' & height: ', rows )
    lefts, rights = [x_min + i * width for i in range(4)], [1 + x_min + i * width for i in range(1, 4)] + [1 + x_min + 4 * width + cols % 4]
    tops, bots = [y_min + i * height for i in range(2)], [y_min + i * height for i in range(1, 2)] + [y_min + 2 * height + rows % 2]

    # Cut the polygon into 8 subtiles
    ROIs = [] ; Tiles = []
    for i in range(4):
        for j in range(2):
            BL = (lefts[i], bots[j])
            TL = (lefts[i], tops[j])
            TR = (rights[i], tops[j])
            BR = (rights[i], bots[j])
            roi = shPol([BL, TL, TR, BR])
            tile = MultiPol.intersection(roi)

            ROIs.append(roi)
            Tiles.append(tile)

    print('Lefts: ', lefts) ; print('Rights: ', rights) ; print('Tops: ', tops) ; print('Bots: ', bots)

    # clean up tiles:
    Tiles_exter = [] ; Tiles_inter = []
    for Ti in Tiles:
        if type(Ti) == shPoint:
            Tiles_exter.append([]) ; Tiles_inter.append([])

        elif type(Ti) == shPol:
            exter = np.asarray(Ti.exterior.xy).T
            Tiles_exter.append([exter])

            if len(Ti.interiors) > 0:
                inters = Ti.interiors
                inter = []
                for i in inters:
                    i = np.asarray(i.xy).T
                    inter.append(i)
                Tiles_inter.append([inter])
            else: Tiles_inter.append([])

        elif type(Ti) == shMPol:
            MPol_list = list(Ti)
            Tiles_exter.append([])
            Tiles_inter.append([])
            for pol in MPol_list:
                exter = np.asarray(pol.exterior.xy).T
                Tiles_exter[-1].append(exter)

                if len(pol.interiors) > 0:
                    inters = pol.interiors
                    inter = []
                    for i in inters:
                        i = np.asarray(i.xy).T
                        inter.append(i)
                    Tiles_inter[-1].append(inter)
                else:
                    Tiles_inter[-1].append([])

        else: Tiles_inter[-1].append([]); Tiles_exter.append(exter)

    # multiprocess to rasterize all tiles
    print('Rasterizing all tiles...')
    pool = mp.Pool(int(mp.cpu_count()-2))
    print('Inner_flag: ' + str(inner_flag))
    print(len(Tiles_exter))
    if inner_flag:
        Tiles_rasterized = pool.starmap(localRasterize, [(Tiles_exter[i], Tiles_inter[i], 1, False, x_min, x_max, y_min, y_max, True) for i in range(8)])
    else: Tiles_rasterized = pool.starmap(localRasterize, [(Tiles_exter[i], [], 1, False, x_min, x_max, y_min, y_max, True) for i in range(8)])
    print('Ready! Merging tiles...')
    pool.close()

    # merge tiles again
    Raster = np.sum(np.asarray(Tiles_rasterized), axis = 0)
    Raster[Raster > 0] = 1

    # print time it took
    tit = time.time() - t_start
    min = int(tit//60)
    sec =  int(tit%60)
    print('It took %d minutes and %d seconds to rasterize the entire polygon!' % (min, sec))

    if get_indices: return(Raster, x_coordinates[:-1], y_coordinates[:-1])
    else: return(Raster)

def bwdist(image):
    """ function to calculate a distance map of a boolean input map (distance to closest cell with value 1)"""
    a = scipy.ndimage.morphology.distance_transform_edt(image==0)
    return(a)

def neigh(image):
    """ neighboorhood function to calculate number of neighboring cells containing value 1
        input has to be boolean"""

    # Neighboorhood alternative
    dim = image.shape
    rows= dim[0]
    cols= dim[1]
    tif_neigh = np.zeros(dim)
    neigh = np.zeros([rows-2,cols-2])
    
    for i in range(0,3):
        for j in range(0,3):
            neigh = neigh+image[i:i+rows-2,j:j+cols-2]
    tif_neigh[1:rows-1,1:cols-1]=neigh 
    tif_neigh = tif_neigh-image
 
    return(tif_neigh)

def rescale(data,rescale_to):
    """ function to rescale linearally"""
    rescaled = rescale_to*(data-np.min(data))/(np.max(data)-np.min(data))
    return(rescaled)

def exportRaster(LU,Example,Path_Name):
    '''
    input is a raster file (.rst, .tif), output a 2D array
    '''
    Example = gdal.Open(Example, gdal.GA_ReadOnly)

    gdal.AllRegister()
    # get info of example data
    rows=Example.RasterYSize
    cols=Example.RasterXSize
    # create output image
    driver=Example.GetDriver()
    outDs=driver.Create(Path_Name,cols,rows,1,GDT_Int32)
    outBand = outDs.GetRasterBand(1)
    outData = np.zeros((rows,cols), np.int16)
    # write the data
    outBand.WriteArray(LU, 0, 0)
    # flush data to disk, set the NoData value and calculate stats
    outBand.FlushCache()
    outBand.SetNoDataValue(-1)

    # georeference the image and set the projection
    outDs.SetGeoTransform(Example.GetGeoTransform())
    outDs.SetProjection(Example.GetProjection())

    del outData

def importRaster(file):
    Raster = gdal.Open(file)
    Band=Raster.GetRasterBand(1)
    Array=Band.ReadAsArray()
    return(Raster, Band, Array)

def coorList2Tiles(XY, rxc = [2,4], add_corner = False):
    """ Method to divide a list of coordinate pairs into a number of tiles to enable parallel processing
    @params:
    XY      -       Required: list of coordinate pairs
    rxc     -       Optional: rows and columns for the tile table
    add_corner      Optional: add the corner of the tiles as coor pair
    """

    pol = shPol(XY)

    r,c = rxc

    r_per_bots = np.arange(0,100,100/r)
    r_per_tops = np.arange(100/r,101,100/r)
    c_per_lefts = np.arange(0,100,100/c)
    c_per_rights = np.arange(100/c,101,100/c)

    XY = np.asarray(XY)
    X,Y = XY[:,0], XY[:,1]
    lefts, rights = np.percentile(X,c_per_lefts), np.percentile(X,c_per_rights)
    bots, tops = np.percentile(Y,r_per_bots), np.percentile(Y,r_per_tops)
    rights[-1] += abs(rights[-1])*0.00001 ; tops[-1] += abs(tops[-1])*0.00001

    Tiles = []
    for i in range(c):
        for j in range(r):
            mask = (X >= lefts[i])*(X<rights[i]) * (Y >= bots[j])*(Y<tops[j])
            XY_tile = XY[mask, :]

            Corners = [(lefts[i], tops[j]),
                        (lefts[i], bots[j]),
                        (rights[i], tops[j]),
                        (rights[i], bots[j])];


            if add_corner:
                Corners_incl = []
                for c in Corners:
                    print(c)
                    if pol.contains(shPoint(c)):
                        Corners_incl.append(c)

                if len(Corners_incl) > 0:
                    C = np.asarray(Corners_incl)
                    print(np.shape(C))
                    print(np.shape(XY_tile))
                    XY_tile = np.vstack((XY_tile, C))



            Tiles.append(XY_tile)

    return Tiles

def image2Blocks(image, vert_n, hor_n, topleftcorners = False):
    """ Function to divide a 2d numpy  array into blocks
    Input: image --> a 2D numpy array
    vert_n--> number of blocks in de vertical
    hor_n --> number of blocks in the horizontal

    Output: Blocks --> list of smaller 2d numpy arrays from the original image
    image_blocks --> numpy 2d array with number of block"""

    rows, cols = np.shape(image)

    Blocks = []
    max_row = 0
    max_col = 0
    t = 1
    TL = []

    image_blocks = np.zeros(np.shape(image))

    for i in range(vert_n + 1):
        for j in range(hor_n + 1):
            height = rows//vert_n
            width = cols//hor_n

            top = i*height
            bottom = (i+1) * height
            left = j*width
            right = (j+1) * width

            block = image[top:bottom,left:right]
            image_blocks[top:bottom,left:right] = t
            Blocks.append(block)

            t += 1

            if bottom > max_row:
                max_row = bottom
            if right > max_col:
                max_col = right

            TL.append([top, left])

    if topleftcorners:
        return Blocks, image_blocks, TL
    else:
        return Blocks, image_blocks

def blocks2image(Blocks, blocks_image):
    """ Function to stitch the blocks back to the original image
    input: Blocks --> the list of blocks (2d numpies)
    blocks_image --> numpy 2d array with numbers corresponding to block number

    output: image --> stitched image """

    image = np.zeros(np.shape(blocks_image))

    for i in range(1,int(np.max(blocks_image))):
        ind = np.asarray(np.where(blocks_image==i))

        top = np.min(ind[0, :])
        bottom = np.max(ind[0, :])
        left = np.min(ind[1, :])
        right = np.max(ind[1, :])
        #print('top: {}, bottom: {}, left: {}, right: {}'.format(top, bottom, left, right))

        image[top:bottom+1,left:right+1] = Blocks[i-1]

    return image

def nestedIndicesList(rows, cols):
    """ Function to get a two column numpy array with all possible indices of a 2D numpy array
    Usage: use to run nested loops in parallel

    Input: rows and cols of your 2d array
    Output: list with 2x1 numpy arrays"""

    ind = np.zeros([rows*cols, 2])
    ind[:,0] = np.tile(np.arange(rows), cols)
    ind[:,1] = np.repeat(np.arange(cols), rows)

    return list(ind)

def retrievePixelValue(geo_coord, image, band):

    """Return floating-point value that corresponds to given point."""

    lat = geo_coord[0]
    lon = geo_coord[1]

    image_np = image.GetRasterBand(band)
    image_np = image_np.ReadAsArray()

    height = np.shape(image_np)[0]
    width = np.shape(image_np)[1]

    # the x and y resolution of the pixels, and the rotation of the
    # raster.
    TL_lon, lon_res, x_rot, TL_lat, lat_rot, lat_res = image.GetGeoTransform()

    # calculate coordinates for each cell

    CoorRasterLon = (TL_lon + 0.5 * lon_res) + np.arange(width) * lon_res
    CoorRasterLat = (TL_lat + 0.5 * lat_res) + np.arange(height) * lat_res

    # get indices from lat and lon
    j = int(round(((lon - (TL_lon - lon_res / 2)) / lon_res)))
    i = int(round(((lat - (TL_lat - lat_res / 2)) / lat_res)))

    return image_np[i,j], i, j

def extent2KML(Image,filename):
    """ Function to plot the extent of a geotiff """
    Band = Image.GetRasterBand(1)
    Array = Band.ReadAsArray()
    rows, cols = np.shape(Array)

    TL_lon, x_res, x_rot, TL_lat, y_rot, y_res = Image.GetGeoTransform()
    res = x_res

    kml = simplekml.Kml()
    pol = kml.newpolygon()

    TL = (TL_lon, TL_lat)
    BL = (TL_lon, TL_lat - rows*res)
    BR = (TL_lon + cols*res, TL_lat - rows*res)
    TR = (TL_lon + cols*res, TL_lat)

    pol.outerboundaryis = [TL,BL,BR,TR,TL]

    kml.save(filename)

def contourKML(Image,cn,filename,minlen = 100):
    """ Function to generate a kml file with contours from a certain raster
     @params:
     the gdal image, a plt.contour output, the filename of the eventual kml file,
     the minimum length of a contour (unit is number of nodes) and color for the eventual
      kml lines (html color codes) """
    c_coor = cn.allsegs[0]

    latlon = Indices2LatLon(Image)

    kml = simplekml.Kml()
    Contour_latlon = []
    Contour_latlon.append([])
    t = 0

    for c_list in c_coor:
        if t%1000 == 0:
            print('Super Process:' +str(round(t/len(c_coor),2)*100)+ ' %')
        t += 1
        f = True
        if len(c_list) > minlen:
            for c in c_list:
                y = c[1]
                x = c[0]

                try:
                    lat = latlon[int(round(y)), int(round(x)), 0]
                    lon = latlon[int(round(y)), int(round(x)), 1]

                    if f == True:
                        Contour_latlon[-1].append((lon, lat))
                        f = False
                    else:
                        lon_old,lat_old = Contour_latlon[-1][-1]
                        dist = ((lon_old-lon)**2+(lat_old-lat)**2)**0.5

                        if dist < 0.0005:
                            Contour_latlon[-1].append((lon, lat))
                        else:
                            Contour_latlon.append([])
                            Contour_latlon[-1].append((lon, lat))
                except Exception as e:
                    exc_type, exc_obj, exc_tb = sys.exc_info()
                    fname = os.path.split(exc_tb.tb_frame.f_code.co_filename)[1]
                    print(exc_type, fname, exc_tb.tb_lineno)

            for C in Contour_latlon:
                l = kml.newlinestring()
                l.name = str(t)
                l.coords = C

    print('Saving file...')
    kml.save(filename)

def calcPolyArea(pol):
    """ Function to calculate the area (units = cells) of a polygon
    input: a numpy array with two cols
    output: the area in number of cells """

    x = pol[:,0]
    y = pol[:,1]
    return 0.5*np.abs(np.dot(x,np.roll(y,1))-np.dot(y,np.roll(x,1)))

def autoCanny(image, sigma=0.33):
    """ Function to use canny without tuning (source: pyimagesearch.com) """
    # compute the median of the single channel pixel intensities
    v = np.median(image)

    # apply automatic Canny edge detection using the computed median
    lower = int(max(0, (1.0 - sigma) * v))
    upper = int(min(255, (1.0 + sigma) * v))
    edged = cv2.Canny(image, lower, upper)

    # return the edged image
    return edged

def shapelyPol2PLTPol(shapelypol, facecolor = 'grey', edgecolor = 'red', alpha = 0.5):
    """ Function to translate a shapely polygon to a plt polygon
    !!! so far only for exterior boundaries
    input: a shapely polygon
    output: a patch to be plotted with plt.add_patch"""

    # get coordinate list
    C = np.rot90(shapelypol.exterior.xy)

    # create a polygon from coordinate list
    pol = Polygon(C, closed=True, facecolor = facecolor, edgecolor = edgecolor, alpha = alpha)  # plt polygon

    return pol

def printProgressBar(iteration, total, prefix = '', suffix = '', decimals = 1, length = 100, fill = '█', printEnd = "\r"):
    """
    Call in a loop to create terminal progress bar
    source: https://stackoverflow.com/questions/3173320/text-progress-bar-in-the-console
    @params:
        iteration   - Required  : current iteration (Int)
        total       - Required  : total iterations (Int)
        prefix      - Optional  : prefix string (Str)
        suffix      - Optional  : suffix string (Str)
        decimals    - Optional  : positive number of decimals in percent complete (Int)
        length      - Optional  : character length of bar (Int)
        fill        - Optional  : bar fill character (Str)
        printEnd    - Optional  : end character (e.g. "\r", "\r\n") (Str)
    """
    percent = ("{0:." + str(decimals) + "f}").format(100 * (iteration / float(total)))
    filledLength = int(length * iteration // total)
    bar = fill * filledLength + '-' * (length - filledLength)
    print('\r%s |%s| %s%% %s' % (prefix, bar, percent, suffix), end = printEnd)
    # Print New Line on Complete
    if iteration == total:
        print()

def shapelyInd2shapelyLatLon(pol, latlon, indlist = False):
    ind = np.asarray(pol.exterior.xy)
    ind_latlon = []

    for i in range(np.shape(ind)[1]):
        col = int(ind[0, i])
        row = int(ind[1, i])

        ind_latlon.append(latlon[row,col,:])

    if indlist:
        return shPol(ind_latlon), ind_latlon
    else:
        return shPol(ind_latlon)

def shapely2KML(pol_ext, path, pol_int = [], name = '', invert = False):
    """ Funtion to translate shapely polygons to kml polygon
    @ params:
        pol_ext - Required: shapely polygon in lat and lon of the exterior boundaries
        path    - Required: DIrectory path or file name (stored in working directory) of the kml file
        pol_int - Optional: List of shapely polygons which describe the interior boundaries
        name    - Optional: String of the name of the polygon
        invert  - Optional: trun on True to change x and y
    """

    kml = simplekml.Kml()
    xy = np.asarray(pol_ext.exterior.xy)
    xy = np.rot90(xy)
    if invert == False:
        xy = list(np.flip(xy))

    xy_int = []
    for i in pol_int:
        xy_i = list(np.flip(np.rot90(np.asarray(i.exterior.xy))))
        xy_int.append(xy_i)

    pol = kml.newpolygon(name = name, outerboundaryis = xy, innerboundaryis = xy_int)
    pol.style.polystyle.color = simplekml.Color.changealphaint(50, simplekml.Color.azure)
    kml.save(path)
    print(f'The Google Earth File is saved in {path}')

def extractPolygon2Snap2FromNDVI(NDVI, LatLon, NDVI_tresh = 0, path = 'Polygon2Snap2.kml'):
    """ Function to extract Polygon2Snap2 from a NDVI image
    @params:
        NDVI        - Required: 2d numpy array of NDVI
        LatLon      - Required: 3d (rowsxcolsx2) numpy array with latitudes and longitudes of the NDVI image
        NDVI_tresh  - Optional: value to distinguish water from non water (default is 0)
        path        - Optional: path to store the resulting kml file (default is 'Polygon2Snap2.kml')
    """

    # Close all open figures
    plt.close('all')

    # treshold value to distinguish water based on NDVI
    t = NDVI_tresh
    Tresh = np.zeros(np.shape(NDVI))

    # apply filter to cancel noise without affecting edges
    blur = cv2.bilateralFilter(NDVI, 3, 5, 5)

    # make a binary image
    (t, binary) = cv2.threshold(src=blur * -1,
                                thresh=t,
                                maxval=255,
                                type=cv2.THRESH_BINARY)

    # convert to proper data type
    binary = np.uint8(binary)

    # contour
    contours = cv2.findContours(image=binary, mode=cv2.RETR_LIST, method=cv2.CHAIN_APPROX_SIMPLE)
    contours = contours[0]
    contours_len = []
    print("Found %d objects." % len(contours))

    # list how long each contour is
    for (i, c) in enumerate(contours):
        contours_len.append(len(c))
        # print("\tSize of contour %d: %d" % (i, len(c)))

    # create empty lists to store the polygons
    patches = []
    sh_Pols = []

    # minimum length of contour to generate polygon
    tresh = 100
    # generate polygons
    t = 0
    print('Generating Polygons...')
    for C in contours:
        t += 1
        printProgressBar(t, len(contours))
        if len(C) > tresh:
            # adapt the dimensions
            C = C[:, 0, :]

            # create a polygon
            pol = Polygon(C, closed=True, facecolor='red', edgecolor='red', alpha=0.05)  # plt polygon
            patches.append(pol)

            shpol = shPol(C)  # shapely polygon
            if type(shpol.buffer(0)) is shPol:
                if shpol.buffer(0).length > tresh:
                    sh_Pols.append(shpol.buffer(0))
            else:
                for s in list(shpol.buffer(0)):
                    if s.length > tresh:
                        sh_Pols.append(s)

    # get the major polygon
    sh_Pols_len = []
    for s in sh_Pols:
        sh_Pols_len.append(s.length)

    C_maj = sh_Pols[np.argmax(sh_Pols_len)]

    # get the interiors
    interiors = []
    t = 0


    print('Getting the interior boundaries...')
    for c in sh_Pols:
        #print("\tSize of contour %d: %d" % (i, c.length))
        t += 1
        printProgressBar(t,len(sh_Pols))
        if c.within(C_maj):
            if c.area < C_maj.area:
                interiors.append(c)
    """
    print('Getting the interior boundaries...')

    def checkWithin(c, C_maj):
        if c.within(C_maj):
            if c.area < C_maj.area:
                return c

    pool = mp.Pool(mp.cpu_count())
    interiors = pool.starmap(checkWithin, [(c, C_maj) for c in sh_Pols])
    """

    interiors = [i for i in interiors if i]
    print(f'Got them! There are {len(interiors)} interior polygons')

    # translate to latlon
    interiors_latlon = []

    for i in interiors:
        pol_latlon = shapelyInd2shapelyLatLon(i, LatLon)
        interiors_latlon.append(pol_latlon)

    pol_maj_latlon = shapelyInd2shapelyLatLon(C_maj, LatLon)

    shapely2KML(pol_maj_latlon, pol_int=interiors_latlon, path=path, name='test')

def indices2LatLon(Image):

    """ Function to calculate a 2D numpy array which hold latitude and longitude in each cell
        input: gdal loaded geotiff
        output: 2D numpy array with lat-lon for each cell"""

    Band=Image.GetRasterBand(1)
    Array=Band.ReadAsArray()
    rows, cols = np.shape(Array)

    TL_lon, x_res, x_rot, TL_lat, y_rot, y_res = Image.GetGeoTransform()
    res = x_res

    latlon = np.zeros([rows,cols,2])
    for i in range(rows):
        printProgressBar(i,rows)
        for j in range(cols):
            lat = TL_lat - i*res
            lon = TL_lon + j*res
            latlon[i,j,0] = lat
            latlon[i,j,1] = lon
    return latlon

def shapely2KML_MultiplePolygons(Pol_list, path, invert = False):
    """ Funtion to translate shapely polygons to kml polygon
      @ params:
          pol_list_element - Required: list of lists where the first element of each list is the outerboundaries of a polygon
          path    - Required: DIrectory path or file name (stored in working directory) of the kml file
          invert  - Optional: trun on True to change x and y
    """

    def savePolList(p, pol, xy, kml, xy_init):
        inxy = np.asarray(pol_list_element[p].exterior.xy)
        inxy = np.rot90(inxy)
        inxy = list(np.flip(inxy))

        xy_init.append(inxy)

        pol = kml.newpolygon(outerboundaryis=xy, innerboundaryis=xy_init)
        pol.style.polystyle.color = simplekml.Color.changealphaint(50, simplekml.Color.azure)

    kml = simplekml.Kml()
    i = 0

    for pol_list_element in Pol_list:

        xy = np.asarray(pol_list_element[0].exterior.xy)
        xy = np.rot90(xy)
        if invert == False:
            xy = list(np.flip(xy))

        if len(pol_list_element) == 1:
            pol = kml.newpolygon(outerboundaryis = xy)
            pol.style.polystyle.color = simplekml.Color.changealphaint(50, simplekml.Color.azure)

        if len(pol_list_element)> 1:
            xy_init = []

            for p in range(1, len(pol_list_element)):
                printProgressBar(p, len(pol_list_element))
                inxy = np.asarray(pol_list_element[p].exterior.xy)
                inxy = np.rot90(inxy)
                inxy = list(np.flip(inxy))

                xy_init.append(inxy)

                pol = kml.newpolygon(outerboundaryis=xy, innerboundaryis=xy_init)
                pol.style.polystyle.color = simplekml.Color.changealphaint(50, simplekml.Color.azure)

    kml.save(path)
    print(f'The Google Earth File is saved in {path}')

def shp2shapely(shp):
    """ Function to load a esri shapefile and transform it to a shapely polygon list
    @ params:
    shp     - Required: directory path to the shapefile"""

    Shapely_List = []

    file = ogr.Open(shp)
    shapefile = file.GetLayer(0)
    # first feature of the shapefile

    n_polygons = shapefile.GetFeatureCount()

    for n in range(n_polygons):
        feature = shapefile.GetFeature(n)
        first = feature.ExportToJson()
        first = json.loads(first)

        shp_geom = shape(first['geometry'])  # or shp_geom = shape(first) with PyShp)

        Shapely_List.append(shp_geom)

    return Shapely_List

def addMiddlePoints(p1, p2, n ):
    """ Function to calculate the middle point between two points
    @ params:
    p1 - Required: set of coordinates of first point
    p2 - Required: set of coordinates of second point
    n  - Required: number of divisions"""

    x_1 = p1[0]
    y_1 = p1[1]

    x_2 = p2[0]
    y_2 = p2[1]

    dx = (x_2 - x_1)/n
    dy = (y_2 - y_1)/n

    ps = [[x_1,y_1]]

    for i in range(1,n):
        x = x_1 + i*dx
        y = y_1 + i*dy

        ps.append([x,y])

    ps.append([x_2,y_2])

    return ps

def addVerticesToPolygon(pol, dist_param, shapely = True):
    """ Function to add vertices on linear segments which are further apart than the defined distance parameter
    @params:
    pol         - Required: shapely polygon
    dist_param  - Required: distance treshold from which extra vertices will be added
    shapely     - Optional (default : true): parameter to indicate whether a shapely polygon is given or just a list with exterior coordinates"""

    if shapely:
        # Exterior
        pol_xy = pol.exterior.coords.xy
        pol_xy = np.asarray(pol_xy)
        pol_xy = np.rot90(pol_xy)
    else:
        pol_xy = np.asarray(pol)

    n = np.shape(pol_xy)[0]
    pol_updates = []

    for i in range(n - 1):

        x_1 = pol_xy[i, 0]
        y_1 = pol_xy[i, 1]

        x_2 = pol_xy[i + 1, 0]
        y_2 = pol_xy[i + 1, 1]

        pol_updates.append((x_1, y_1))

        dist = np.sqrt((x_1 - x_2) ** 2 + (y_1 - y_2) ** 2)

        if dist > dist_param:
            n_inter = int(dist // (dist_param / 2))
            inter = addMiddlePoints((x_1, y_1), (x_2, y_2), n_inter)
            for i in inter:
                pol_updates.append(i)
        else:
            pol_updates.append((x_2, y_2))

    pol_updates.append((pol_xy[-1, 0], pol_xy[-1, 1]))
    exterior = pol_updates

    if shapely:
        # Interiors
        interiors = []

        for interior in pol.interiors:
            pol_xy = interior.coords.xy
            pol_xy = np.asarray(pol_xy)
            pol_xy = np.rot90(pol_xy)

            n = np.shape(pol_xy)[0]
            pol_updates = []

            for i in range(n-1):

                x_1 = pol_xy[i, 0]
                y_1 = pol_xy[i, 1]

                x_2 = pol_xy[i+1, 0]
                y_2 = pol_xy[i+1, 1]

                pol_updates.append((x_1, y_1))

                dist = np.sqrt((x_1-x_2)**2+(y_1-y_2)**2)

                if dist > dist_param:
                    n_inter = int(dist//(dist_param/2))
                    inter = addMiddlePoints((x_1,y_1),(x_2,y_2), n_inter)
                    for i in inter:
                        pol_updates.append(i)
                else:
                    pol_updates.append((x_2, y_2))

            pol_updates.append((pol_xy[-1, 0], pol_xy[-1, 1]))
            interiors.append(pol_updates)


        pol_updated = shPol(shell = exterior, holes = interiors)

        return pol_updated
    else:
        return exterior

def snapPolygonstoPolygon(Polygon2Snap, Polygon2Snap2, Snap_Sensitivity, n_blocks = 10, Progressbar = False):
    """ Function to snap a polygon to another polygon
    @ params:
    polygon2Snap        - Required  : Polygon which has to be edited to snap to a reference polygon
    polygon2Snap2       - Required  : Reference Polygon
    Snap_Sensitivity    - Required  : Max distance to snap
    n_blocks            - Optional  : Number of blocks per row """

    Polygon2Snap2_mult = []
    Polygon2Snap2_xy = []

    t = 0

    # divide the polygon to snap to in blocks
    Polygon2Snap2_mult = []
    Polygon2Snap2_xy = []
    Polygon2Snap2_x = []
    Polygon2Snap2_y = []

    for c in Polygon2Snap2:
        if type(c) is shPol:
            Polygon2Snap2_mult.append([c])
            c = addVerticesToPolygon(c, Snap_Sensitivity / 4)
            xy = c.exterior.xy
            xy = np.rot90(xy)
            for i in xy:
                Polygon2Snap2_xy.append(i)
                Polygon2Snap2_x.append(i[0])
                Polygon2Snap2_y.append(i[1])
            for int_pol in c.interiors:
                xy = int_pol.xy
                xy = np.rot90(xy)
                for i in xy:
                    Polygon2Snap2_xy.append(i)
                    Polygon2Snap2_x.append(i[0])
                    Polygon2Snap2_y.append(i[1])

        else:
            for s in list(c):
                Polygon2Snap2_mult.append([s])
                s = addVerticesToPolygon(s, Snap_Sensitivity / 4)
                xy = s.exterior.xy
                xy = np.rot90(xy)
                for i in xy:
                    Polygon2Snap2_xy.append(i)
                    Polygon2Snap2_x.append(i[0])
                    Polygon2Snap2_y.append(i[1])
                ints = s.interiors
                for int_pol in ints:
                    xy = int_pol.xy
                    xy = np.rot90(xy)
                    for i in xy:
                        Polygon2Snap2_xy.append(i)
                        Polygon2Snap2_x.append(i[0])
                        Polygon2Snap2_y.append(i[1])

    Polygon2Snap2_x = np.asarray(Polygon2Snap2_x)
    Polygon2Snap2_y = np.asarray(Polygon2Snap2_y)
    Polygon2Snap2_xy = np.asarray(Polygon2Snap2_xy)

    left = np.min(Polygon2Snap2_x)
    right = np.max(Polygon2Snap2_x)

    bottom = np.min(Polygon2Snap2_y)
    top = np.max(Polygon2Snap2_y)

    dx = (right - left) / n_blocks
    dy = (top - bottom) / n_blocks

    if dx < Snap_Sensitivity:
        print('The number of blocks is too high!')
    elif dy < Snap_Sensitivity:
        print('The number of blocks is too high!')

    Polygon2Snap2_Blocks = []
    Blocks_lefts = np.zeros(n_blocks)
    Blocks_tops = np.zeros(n_blocks)

    if Progressbar: print('Creating the blocks...')
    for j in range(n_blocks):
        for i in range(n_blocks):
            if Progressbar: printProgressBar(j, n_blocks)
            Blocks_lefts[i] = left + i * dx
            Blocks_tops[j] = top - j * dy
            block = Polygon2Snap2_xy[
                    (Polygon2Snap2_x > left + i * dx) & (Polygon2Snap2_x <= left + (i + 1) * dx) & (Polygon2Snap2_y < top - j * dy) & (
                                Polygon2Snap2_y >= top - (j + 1) * dy), :]
            Polygon2Snap2_Blocks.append(block)

    Blocks_lefts = np.asarray(Blocks_lefts)
    Blocks_tops = np.asarray(Blocks_tops)

    # prepare the snapping
    Polygon2Snap = addVerticesToPolygon(Polygon2Snap, Snap_Sensitivity/4)
    Man_snapped = Polygon2Snap
    Man_xy = Polygon2Snap.exterior.xy
    Man_snapped_xy = Man_snapped.exterior.xy
    Man_xy = np.rot90(Man_xy)

    Man_snapped_xy_tuplelist = []

    # carry out the snapping

    if Progressbar: print('Snapping...')
    for i in range(np.shape(Man_xy)[0]):
        if Progressbar: printProgressBar(i,np.shape(Man_xy)[0])

        x = Man_xy[i][0]
        y = Man_xy[i][1]

        # in what block is the point2snap located?
        diff_x = Blocks_lefts - x
        xn = np.where(diff_x < 0, diff_x, -np.inf).argmax()
        diff_y = Blocks_tops - y
        yn = np.where(diff_y > 0, diff_y, np.inf).argmin()
        n = n_blocks * yn + xn

        Sub2SnapTo = []
        for nn in [-(n_blocks+1),-(n_blocks),-(n_blocks-1),-1,0,1,n_blocks-1,n_blocks,n_blocks+1]:
            for xy in Polygon2Snap2_Blocks[n + nn]:
                Sub2SnapTo.append(xy)

        Dist = []

        for j in range(len(Sub2SnapTo)):
            x_ref = Sub2SnapTo[j][0]
            y_ref = Sub2SnapTo[j][1]

            d = np.sqrt((x-x_ref)**2+(y-y_ref)**2)

            Dist.append(d)

        if len(Dist)>0:
            minDist = np.min(Dist)
            minDist_arg = np.argmin(Dist)

            if minDist < Snap_Sensitivity:
                x_snap = Sub2SnapTo[minDist_arg][0]
                y_snap = Sub2SnapTo[minDist_arg][1]
                Man_snapped_xy_tuplelist.append((x_snap,y_snap))
                #print(f'Point Snapped! from x: {round(x,4)} y: {round(y,4)} to from x: {round(x_snap,4)} y: {round(y_snap,4)}')

            else:
                Man_snapped_xy_tuplelist.append((x,y))
        else:
            Man_snapped_xy_tuplelist.append((x, y))


    Man_snapped_xy_tuplelist.append(Man_snapped_xy_tuplelist[0])
    Man_snapped = shPol(Man_snapped_xy_tuplelist)

    return Man_snapped

def KML2Grid(Input_fn,Output_fn, Points_resolution = None):
    """ Function to transform KML file into a meshed grid
    @ params:
    Input_fn            - Required  : Directory to indicate input google earth kml file
    Output              - Required  : Directory to indicate name and location output vtk file (readable by gmsh software)
    Points_Resolution   - Optional  : List of lists where each list contains Lat Lon and resolution """

    # read kml as multi string
    kml_str = open(Input_fn, 'r', encoding='latin-1').read()
    kml_splitted = np.asarray(kml_str.splitlines())
    loc = np.flatnonzero(np.core.defchararray.find(kml_splitted, '<coordinates>') != -1)

    # parse lat and lon into lists
    t = 0
    coor_int = []
    for l in loc:
        if t == 0:
            coor = kml_splitted[l + 1]
            coor = np.asarray(re.split(' |/t|,', coor))

        else:
            c_int = kml_splitted[l + 1]
            c_int = np.asarray(re.split(' |/t|,', c_int))
            coor_int.append(c_int)

        t += 1

    Coor_UTM = []
    for i in range(0, len(coor) - 1, 3):
        try:
            x, y, utm_number, utm_letter = utm.from_latlon(float(coor[i + 1]), float(coor[i]))
            Coor_UTM.append([x, y, 0])
        except:
            pass

    Coor_UTM_int = []
    for c_int in coor_int:
        C_UTM_int = []
        for i in range(0, len(c_int) - 1, 3):
            try:
                x, y, utm_number, utm_letter = utm.from_latlon(float(c_int[i + 1]), float(c_int[i]))
                C_UTM_int.append([x, y, 0])
            except:
                pass
        Coor_UTM_int.append(C_UTM_int)

    ext = Coor_UTM
    int = Coor_UTM_int

    # create polygon
    geom = pygmsh.built_in.Geometry()

    geom_points = []

    for p in ext:
        geom_points.append(geom.add_point(p, lcar=10000))

    geom_points.append(geom_points[0])
    geom_spline = geom.add_bspline(geom_points)
    geom_ll = geom.add_line_loop([geom_spline])
    geom_surf = geom.add_plane_surface(geom_ll)

    mesh = pygmsh.generate_mesh(geom, dim=2)

    # store mesh by means of meshio

    fn = Output_fn
    if fn[-4:] == '.vtk':
        pass
    else:
        fn = fn + '.vtk'
    meshio.write(fn, mesh)

    print(fn + ' is saved!')

    return ext,int

def polPlot(XY, XY_inner = None, plottitle = None, showvertices = False, showverticeslabels = False, showverticeslabelsinterval = 1, plotonaxis = 0, empty = False, vertices_color = 'silver'):
    """ Function to plot polygon based on list of coordinates
    XY      -       Required : List of x and y coordinates (list of lists)
    """
    if plotonaxis == 0:
        f, a = plt.subplots()
    else:
        a = plotonaxis

    XY_np = np.asarray(XY)
    xmin = np.min(XY_np[:, 0])
    xmax = np.max(XY_np[:, 0])
    ymin = np.min(XY_np[:, 1])
    ymax = np.max(XY_np[:, 1])

    if plottitle:
        a.set_title(plottitle, fontweight = 'bold')

    a.set_xlim(xmin - 0.1*(xmax-xmin), xmax + 0.1*(xmax-xmin))
    a.set_ylim(ymin - 0.1*(ymax-ymin), ymax + 0.1*(ymax-ymin))

    if empty:
        col = 'none'
        colin = 'none'
    else:
        col = 'silver'
        colin = 'white'

    pol = Polygon(XY, facecolor = col, edgecolor = 'turquoise' ) # matplotlib.patches.Polygon

    if XY_inner is not None:
        pols_inner = []
        for i in XY_inner:
            pol_inner = Polygon(i, facecolor=colin, edgecolor='turquoise')  # matplotlib.patches.Polygon
            pols_inner.append(pol_inner)

    a.add_patch(pol)
    if XY_inner is not None:
        for pol_inner in pols_inner:
            a.add_patch(pol_inner)

    if showvertices:
        a.scatter(XY_np[:,0], XY_np[:,1], s = 4, edgecolor = vertices_color, facecolor = 'none', zorder = 10)
        if XY_inner:
            for i in XY_inner:
                a.scatter(np.asarray(i)[:, 0], np.asarray(i)[:, 1], s=4, edgecolor= vertices_color, facecolor='none', zorder = 10)

        t = 1

        if showverticeslabels:
            for x,y in XY:
                if t%showverticeslabelsinterval == 0: a.annotate(str(t),(x,y))
                t += 1
            if XY_inner:
                for i in XY_inner:
                    for x,y in i:
                        if t%showverticeslabelsinterval == 0: a.annotate(str(t), (x, y))
                        t += 1


    a.set_aspect('equal')

def loadPolygonFromShapefile(fn, Print_Coordinate_System = False):
    """ Function to load a shapefile containing one polygon with outer and inner boundaries
    fn      -       Required: File path directory
    """

    # load shapefile with the ogr toolbox of osgeo
    file = ogr.Open(fn)
    shape = file.GetLayer(0)
    if Print_Coordinate_System:
        proj = str(shape.GetSpatialRef())
        proj = proj.split(',', 2)[0][8:-1]
        print(f'The Coordinate system is: {proj}')

    n_features = shape.GetFeatureCount()

    assert n_features == 1

    for i in range(n_features):
        feature = shape.GetFeature(i)
        feature_JSON = feature.ExportToJson()
        feature_JSON = json.loads(feature_JSON)

        coor = feature_JSON['geometry']['coordinates']

    if len(coor[0][0]) > 3:
        coor = coor[0]
    exterior = coor[0]
    interiors = coor[1:]

    return exterior, interiors, feature_JSON

def simplifyPol(XY,tresh):
    """ Function to remove nodes which are located nearby other nodes
    XY      -       Required: coordinate list of nodes
    tresh   -       Required: treshold distance
    """


    def remove_below_tresh(XY,tresh):

        n = len(XY)

        def eucldist(p1,p2):
            x1 = p1[0]
            x2 = p2[0]
            y1 = p1[1]
            y2 = p2[1]

            if len(p1) >2:
                z1 = p1[2]
                z2 = p2[2]
                ed = ((x2-x1)**2 + (y2-y1)**2 + (z2 -z1)**2)**0.5
            else:
                ed = ((x2 - x1) ** 2 + (y2 - y1) ** 2) ** 0.5
            return ed

        d_l = eucldist(XY[-1],XY[0])
        dist_l = []
        dist_l.append(d_l)

        d_r = eucldist(XY[0],XY[1])
        dist_r = []
        dist_r.append(d_r)

        for i in range(1,n-1):
            d_l = eucldist(XY[i-1], XY[i])
            dist_l.append(d_l)

            d_r = eucldist(XY[i], XY[i+1])
            dist_r.append(d_r)

        d_l = eucldist(XY[-2], XY[-1])
        dist_l.append(d_l)

        d_r = eucldist(XY[-1], XY[0])
        dist_r.append(d_r)

        el = np.zeros([n,1])
        for i in range(n):
            if dist_l[i] < tresh and dist_r[i] < tresh:
                el[i] = 1

        XY_fil = []
        odd = True
        for i in range(n):
            if el[i] == 0:
                XY_fil.append(XY[i])
            elif el[i] == 1 and odd:
                odd = False
            elif el[i] == 1 and odd == False:
                odd = True
                XY_fil.append(XY[i])
            else:
                pass

        return XY_fil

    XY_n = remove_below_tresh(XY, tresh)

    while len(XY_n) != len(XY) and len(XY) > 4:
        XY = XY_n
        XY_n = remove_below_tresh(XY, tresh)

    return XY_n

def WGS842UTM(XY, utm_n, utm_l):
    """Method to transform a list of lat lon coordinates in WGS84 to a utm system
    XY      -   Required: list of coordinate pairs
    utm_n   -   Required: UTM zone number
    utm_l   -   Required: UTM zone letter """

    XY_utm = []

    for lonlat in XY:
        lon, lat = lonlat[0:2]
        utmxy = utm.from_latlon(lat, lon, utm_n, utm_l)
        XY_utm.append(utmxy[0:2])

    return XY_utm

def removeIsolations(XY, tresh_dist, plot_progress=False):
    """ Method to remove isolated points/vertices from a list
    @params:
    XY          -   Required: list of coordinate pairs
    tresh_dist  -   Required: distance to distinguish isolated cells from non-isolated cells
    """

    XY_np = np.asarray(XY)
    XY_upd = []

    t = 0
    for x, y in XY:

        t += 1
        if plot_progress:
            printProgressBar(t, len(XY))
        x_l = x - tresh_dist
        x_r = x + tresh_dist
        y_b = y - tresh_dist
        y_t = y + tresh_dist

        mask = (XY_np[:, 0] > x_l) * (XY_np[:, 0] < x_r) * (XY_np[:, 1] > y_b) * (XY_np[:, 1] < y_t)
        XY_np_close = XY_np[mask, :]

        XY_neigh = []
        for i in range(len(XY_np_close[:, 0])):
            x_ = XY_np_close[i, 0]
            y_ = XY_np_close[i, 1]
            dist = np.sqrt((x - x_) ** 2 + (y - y_) ** 2)
            if dist <= tresh_dist:
                XY_neigh.append([x_, y_])

        if len(XY_neigh) > 1:
            XY_upd.append([x, y])

    return XY_upd

def makeNeighbouringLists(XY, tresh_dist, plot_progress=False):
        """ Method to get neighbouring lists
        @params:
        XY          -   Required: list of coordinate pairs
        tresh_dist  -   Required: distance to distinguish isolated cells from non-isolated cells
        """

        XY_np = np.asarray(XY)
        XY_upd = []
        XY_closest = []
        XY_dist = []

        t = 0
        for x, y in XY:

            t += 1
            if plot_progress:
                printProgressBar(t, len(XY))
            x_l = x - tresh_dist
            x_r = x + tresh_dist
            y_b = y - tresh_dist
            y_t = y + tresh_dist

            mask = (XY_np[:, 0] > x_l) * (XY_np[:, 0] < x_r) * (XY_np[:, 1] > y_b) * (XY_np[:, 1] < y_t)
            XY_np_close = XY_np[mask, :]

            XY_neigh = []
            XY_neigh_dists = []

            for i in range(len(XY_np_close[:, 0])):
                x_ = XY_np_close[i, 0]
                y_ = XY_np_close[i, 1]
                dist = np.sqrt((x - x_) ** 2 + (y - y_) ** 2)
                if dist <= tresh_dist and dist > 0:
                    XY_neigh.append([x_, y_])
                    XY_neigh_dists.append(dist)

            if len(XY_neigh) == 0:
                XY_closest.append([x, y])
                XY_dist.append(0)
            else:
                xy_closest = XY_neigh
                xy_closest_dist = XY_neigh_dists
                XY_closest.append(xy_closest)
                XY_dist.append(xy_closest_dist)

        return XY_closest, XY_dist

def localMerge(treshold, XY):
    """  local function to feed starmap
    """
    XY = list(XY)
    # copy original method
    XY_orig = XY.copy()
    # keep length of original list
    Len_orig = len(XY)
    # make sure the while loop will start
    Len_pre = Len_orig + 1

    # counter
    t = 0

    # keep removing untill there is nothing more to remove
    while Len_pre != len(XY):
        # save the length of the coordinate list before possibly removing one pair
        Len_pre = len(XY)
        # create a list per point with all neighbours and their corresponding distances
        XY_closest, XY_dist = makeNeighbouringLists(XY, treshold)

        # convert to numpy arrays
        XY_np = np.asarray(XY)
        XY_closest_np = np.asarray(XY_closest)

        # create a list for each point with the closest neighbours (both the index and the distance)
        closest_dist = []
        closest_argdist = []
        for xy in XY_dist:
            closest_dist.append(np.min(xy))
            closest_argdist.append(np.argmin(xy))
        # convert to numpy array
        closest_dist = np.asarray(closest_dist)
        # vertices without number should get a very large distance so they would never be the minimum
        closest_dist[closest_dist == 0] = np.max(closest_dist) * 1000

        # get closest neighbour and its index
        closest_min = np.min(closest_dist)
        closest_argmin = np.argmin(closest_dist)

        # only proceed if there is a neighbour selected
        if closest_min < np.max(closest_dist) * 1000:
            P1 = XY[closest_argmin]
            P2 = XY_closest[closest_argmin][closest_argdist[closest_argmin]]
            P_m = np.mean(np.vstack((P1, P2)), axis=0)

            XY_upd = np.asarray(XY)
            X,Y = XY_upd[:,0], XY_upd[:,1]

            P1x, P1y = P1
            P2x, P2y = P2

            X_mask = (X == P1x) | (X == P2x)
            Y_mask = (Y == P1y) | (Y == P2y)
            mask = np.invert(X_mask)*np.invert(Y_mask)
            XY_upd = list(XY_upd[mask,:])
            XY_upd.append(list(P_m))

            XY = XY_upd
        t += 2

    return (XY)

def mergeCloseNeighbours(XY, tresh, print_orig = False):
    """ Method to replace close neighbouring points with one
    @params:
    XY          -   Required: list of coordinate pairs
    tresh       -   Required: distance to distinguish isolated cells from non-isolated cells
    """

    XY_orig = XY.copy()
    # create tiles
    Tiles = coorList2Tiles(XY, rxc=[2, 3])

    # starmap local function
    pool = mp.Pool(mp.cpu_count())
    #Tiles_Neigh = pool.map(functools.partial(localMerge, tresh), Tiles)
    Tiles_Neigh = pool.starmap(localMerge, [(tresh, xy) for xy in Tiles])
    pool.close()

    xy_np = []
    for t in Tiles_Neigh:
        xy_np.append((np.asarray(t)))
    res = np.vstack(xy_np)

    res = localMerge(tresh, res)

    if print_orig:
        return res, XY_orig
    else:
        return res

def getSkeletonPoints(XY_exter, Iso_tresh, XY_inter = []):
    """ Method to calculate the skeleton points of a polygon
    @params:
    XY_exter    -   Required: List of points which form the outer boundary of a polygon
    Iso_tresh   -   Required: Treshold distance to determine isolated points
    XY_inter    -   Optional: List of Lists of points which form the inner boundaries of a polygon
    """

    # create coordinate list with all coordinates
    XY_all = XY_exter.copy()
    if len(XY_inter)>0:
        for i in XY_inter:
            XY_all += i

    # create Voronoi points
    vor = Voronoi(XY_all)
    """voronoi_plot_2d(vor)"""
    vor_vert = vor.vertices

    XY_inter_upd = []
    for xy in XY_inter:
        if len(xy) >= 3:
            XY_inter_upd.append(xy)
    XY_inter = XY_inter_upd

    # create shapely polygon
    if len(XY_inter)==0:
        pol = shPol(XY_exter)
    else:
        pol = shPol(XY_exter, XY_inter)
    Vor_int = []

    # remove all points located outside the polygon
    for i in range(len(vor_vert)):
        #printProgressBar(i, len(vor_vert))
        p = shPoint(vor_vert[i, 0], vor_vert[i, 1])
        if p.within(pol):
            Vor_int.append([vor_vert[i, 0], vor_vert[i, 1]])

    # remove isolated points and merge close neighbours
    Vor_int_upd = removeIsolations(Vor_int, Iso_tresh, plot_progress=True)
    Vor_int_upd, Vor_int_upd_orig = mergeCloseNeighbours(Vor_int_upd, Iso_tresh, print_orig=True)
    Vor_int_upd2 = removeIsolations(Vor_int_upd, Iso_tresh * 2, plot_progress=True)
    Vor_int_upd3, Vor_int_upd2 = mergeCloseNeighbours(Vor_int_upd2, Iso_tresh * 2, print_orig=True)

    return Vor_int_upd3, pol

def createMergeResolutionField(xy, outputfile, minres = 25, zeropoint = [0,0], xres = 1, yres = 1):
    """
        TELEMAC: Function to write a field file to accompagny a .geo file to determine resolution

    :param xy: 2d numpy array with the values of mesh size
    :param outputfile: directory to store the newly generated file
    :param minres: minimum mesh size
    :param zeropoint: coordinate pair of the left bottom corner of the numpy image
    :param xres: cell width
    :param yres: cell height
    :return: a .bg file
    """

    print('Writing Mesh Size Field File...')

    xy[xy < minres] = minres

    f = open('Scratch/header_file.txt', "w+")
    f.write('%f %f %f\n' % (zeropoint[0] + xres / 2, zeropoint[1] + yres / 2, 0))
    f.write('%f %f %f\n' % (xres, yres, 1))
    f.write('%d %d %d\n' % (np.shape(xy)[1], np.shape(xy)[0], 1))
    f.close()

    xy = np.flip(xy.T, axis=1).flatten()
    n = len(xy)//8

    def writeFile(n, lis):
        f = open('Scratch/file_%s.txt' % (str(n)), "w+")
        t = 0
        for v in lis:
            t += 1;
            f.write('%f\n' % (v))
        f.close()
        return 'Scratch/file_%s.txt' % (str(n))

    n, n_rest = len(xy) // 8, len(xy) % 8
    N, N[0], N[-1] = [n for i in range(9)] * np.arange(9), 0, len(xy)
    xy_tiles = [xy[N[i]:N[i + 1]] for i in range(8)]

    pool = PathosPool(8)
    P = pool.map(writeFile, np.arange(8), xy_tiles)
    P.insert(0,'Scratch/header_file.txt')


    with open(outputfile, 'wb') as wfd:
        for f in P:
            with open(f, 'rb') as fd:
                shutil.copyfileobj(fd, wfd)

    for f in P: os.remove(f)

    print('Mesh Size Field File is written!')

def arr2Geotiff(arr, fname, TL, Res, projEPSG):
    """
    Function to export a 2D numpy array as geotiff
    :param arr: 2D numpy array
    :param fname: output file directory for tiff file
    :param TL: tuple or list of x and y coordinate of the top left corner
    :param Res: tuple, list of single number with resolution
    :param projEPSG: EPSG code
    :return:
    """

    # Transform array to 64 bit floats (readable in QGIS)
    arr = np.float64(arr)
    # get array dimensions
    rows, cols = np.shape(arr)

    # Create gdal driver in the geo-tiff format
    drv = gdal.GetDriverByName("GTiff")
    ds = drv.Create(fname, cols, rows, 1, gdal.GDT_Float64)

    # get coordinates of topleft corner
    TL_x, TL_y = TL

    # get resolution
    if type(Res) == list:
        x_res, y_res = Res
    elif type(Res) == tuple:
        x_res, y_res = Res
    else:
        x_res, y_res = Res, Res

    # set geometry data (TL coordinates, resolution and skew)
    args = (TL_x, x_res, 0, TL_y, 0, y_res)
    ds.SetGeoTransform(args)

    # write band 1 values
    ds.GetRasterBand(1).WriteArray(arr)
    ds.GetRasterBand(1).SetNoDataValue(-9999)

    # set projection based on given EPSG
    srs = osr.SpatialReference()
    srs.ImportFromEPSG(projEPSG)
    ds.SetProjection(srs.ExportToWkt())

    # Close
    ds = None

def unique(list1):
    """
    Function to get unique values of a list in same order
    :param list1: list
    :return:  list
    """

    # make sure it is a list
    list1 = list(list1)

    # intilize a null list
    unique_list = []

    # traverse for all elements
    for x in list1:
        # check if exists in unique_list or not
        if x not in unique_list:
            unique_list.append(x)

    return unique_list

def removeVertices(xy, geotype = 'pol', sens = 0):
    """ Function to remove excessive vertices of line/polygon
    :param xy: vertices list of a polygon (list of pairs or nx2 numpy array with x- and y-coordinates of the polygon's vertices)
    :param geotype: 'pol' or 'line'
    :param sens: maximum allowed change of slope
    :return: vertices list without excessive points (nx2 numpy array with x- and y-coordinates of the polygon's vertices)
    """

    # transform list to np array
    if type(xy) == list: xy = np.asarray(xy)

    # get number of vertices
    n = np.shape(xy)[0]

    # initialize an empty numpy array to store the items to be removed
    rem = np.zeros(n)

    # In case of a line, never remove the first or last vertex, in case of polygon that is possible
    if geotype ==  'pol': start_loop, end_loop = 0, n
    else: start_loop, end_loop = 1, n-1

    # loop over all vertices
    for i in range(start_loop, end_loop):
        # get left and right neighbor indices
        if i == 0: l, r = -1,1
        elif i == n-1: l, r = n-2, 0
        else: l, r = i-1, i+1

        # check slope of straights
        x, y = xy[i,:]
        lx, ly = xy[l,:]
        rx, ry = xy[r,:]
        la = (ly - y)/(lx - x)
        ra = (ry - y)/(rx - x)

        # if slopes are equal, the vertex is not essential
        if abs(la - ra) <= sens:
            rem[i] = 1

    # mask to remove all excessive vertices
    xy_upd = xy[rem==0]
    return xy_upd









