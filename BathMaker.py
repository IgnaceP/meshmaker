import numpy as np
import scipy
import matplotlib.pyplot as plt
from scipy.interpolate import griddata
import pandas as pd
from shapely.geometry import Polygon as shPol
from shp2mpol import *
from pol2pol import *

# Parameter Definition
#############################################################################################################
buff_size = 5
elev = 2.5

# prescribed angle for the slope of the banks (tan)
angle = 0.001

# export the bathymetric map as geotiff
tif_flag = True
tif_fn = 'Scratch/test.tif'

# plot resulting scatter distribution?
plot_flag = False

channels = '/home/ignace/Documents/PhD/Data/Land Use/Guayas_Land_Use/Churute_Subregion_Channels.shp'
basin = '/home/ignace/Documents/PhD/Data/Land Use/Guayas_Land_Use/Churute_Subregion.shp'
outputfile = '/home/ignace/Documents/PhD/TELEMAC/Playing_Around/First_Steps/Test_Cases/input/v3/Churute_Subregion.bath'
#############################################################################################################

# Load Polygon
pol, epsg = shp2Mpol(channels, return_coordinate_system = True)
pol = pol2Pol(pol, orig_epsg= epsg, new_epsg= 32717)[0]

# convert to Numpy Array and rotate 90° do the shape becoms n x 2, remove z-coordinates
exter = np.rot90(pol.exterior.xy)[:-1,0:2]
inter = [np.rot90(i.xy)[:-1,0:2] for i in pol.interiors]

# create buffer for stable slope
pol_buff = pol.buffer(buff_size)
buff_exter = np.rot90(np.asarray(pol_buff.exterior.xy))
buff_inners = []
for i in range(len(inter)):
    buff_inners.append(np.rot90(np.asarray(pol_buff.interiors[i].xy)))
# fun.polPlot(buff_exter, XY_inner=buff_inners)

# create 2d numpy arrays with all coordinate pairs from the buffer
buff_xy = buff_exter
for buff_i in buff_inners:
    buff_xy = np.vstack((buff_xy,buff_i))

buff_z = np.zeros([np.shape(buff_xy)[0],1]) + elev
buff_xy = np.hstack((buff_xy, buff_z))

# create 2d numpy arrays with all coordinate pairs from the channels
chan_xy = exter
chan_xy_pol = exter
for i in inter:
    chan_xy = np.vstack((chan_xy, i))
chan_z = np.zeros([np.shape(chan_xy)[0], 1])
chan_xy = np.hstack((chan_xy, chan_z))

# create 2d numpy arrays with all coordinate pairs from the corners
exter, inter, json = fun.loadPolygonFromShapefile(basin)
exter = np.asarray(fun.WGS842UTM(exter, 17, 'M'))

x_min = np.min(exter[:,0])-1000
x_max = np.max(exter[:,0])+1000
y_min = np.min(exter[:,1])-1000
y_max = np.max(exter[:,1])+1000

corners = np.zeros([4,3])
corners[0,:] = np.asarray([x_min,y_min,elev])
corners[1,:] = np.asarray([x_min,y_max,elev])
corners[2,:] = np.asarray([x_max,y_min,elev])
corners[3,:] = np.asarray([x_max,y_max,elev])

# Add Distancemap
try:
    ras, X_orig, Y_orig = np.load('Scratch/ras.npy'),np.load('Scratch/X_orig.npy'), np.load('Scratch/Y_orig.npy')
except:
    ras, X_orig, Y_orig = fun.localRasterizeInParallel(chan_xy_pol, XY_inner_L=[], generalize_factor= 10, get_indices= True)
    np.save('Scratch/ras.npy', ras)
    np.save('Scratch/X_orig.npy', X_orig)
    np.save('Scratch/Y_orig.npy', Y_orig)

x, y = np.meshgrid(X_orig, Y_orig); x, y = x.flatten(), y.flatten()
dist = (fun.bwdist(ras)*10)

buff = np.zeros(np.shape(dist))
buff[dist<= 5] = 1

top = angle*np.max(dist)
dist = (dist - np.min(dist))/np.max(dist)*top + elev

if tif_flag:
    TL = np.min(X_orig), np.max(Y_orig)
    res = Y_orig[1] - Y_orig[0]
    arr2Geotiff(dist, tif_fn, TL, res, 32317)
    print('Geotiff file created as %s.' % (tif_fn))

dist_fl = dist.flatten()

x_dist = np.asarray([x[i] for i in range(len(x)) if dist_fl[i] > elev and i%25 == 0])
y_dist = np.asarray([y[-i] for i in range(len(x)) if dist_fl[i] > elev and i%25 == 0])
z_dist = np.asarray([dist_fl[i] for i in range(len(x)) if dist_fl[i] > elev and i%25 == 0])

# stack all arrays
XYZ = np.vstack((chan_xy, buff_xy, corners))
X = np.hstack((XYZ[:,0],x_dist))
Y = np.hstack((XYZ[:,1],y_dist))
Z = np.hstack((XYZ[:,2],z_dist))

if plot_flag:
    plt.scatter(X, Y, 10, Z)

df = pd.DataFrame({'X':X, 'Y':Y, 'Z':Z})
df.to_pickle(outputfile)
print('Bath file created as %s.' % (outputfile))

print('Ready!')
#--------------------------------------------------------------------------
