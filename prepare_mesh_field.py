""" Python script to generate a input file for a structured field in GMSH based on a 2d numpy array """

import numpy as np
import matplotlib.pyplot as plt
import Scripts.Functions as fun
from scipy.interpolate import griddata
import pickle

###***** Churute *****###

# Channels
#############################################################################################################
x_orig = np.load('ChannelWidth_Churute/rasterized_Channels_Churute_25m_Xorig.npy')
y_orig = np.load('ChannelWidth_Churute/rasterized_Channels_Churute_25m_Yorig.npy')
bandwidth = 250

# get zeropoint and resolution
zeropoint = [min(x_orig)-bandwidth*25,min(y_orig)-bandwidth*25]
xres = x_orig[1]-x_orig[0]
yres = y_orig[1]-y_orig[0]

# load channel width
xy = np.load('ChannelWidth_Churute/ChannelWidth_Churute_25m_v1.npy')
xy = np.pad(xy, pad_width=bandwidth, mode='constant', constant_values=-100)
rows, cols = np.shape(xy)


xi, yi = np.arange(0,cols), np.arange(0, rows)
xi, yi = np.meshgrid(xi,yi)
X, Y, D = xi.flatten(), yi.flatten(), xy.flatten()
X, Y, D = X[D>=0], Y[D>=0], D[D>=0]
grid = griddata((X,Y), D, (xi, yi), method = 'nearest')

mangroves = (xy == - 100); channels = (xy >= 0)
mangroves_dist = fun.bwdist(channels)*25

xy[mangroves] = mangroves_dist[mangroves] + grid[mangroves]
xy = xy/4

pathname = '/home/ignace/Documents/PhD/TELEMAC/Playing_Around/GMSH/Churute_MeshSize_NumpyArray.pck'

XY = [xy, zeropoint, xres, yres]
with open(pathname, 'wb') as output:
    pickle.dump(XY, output, pickle.HIGHEST_PROTOCOL)


#############################################################################################################
