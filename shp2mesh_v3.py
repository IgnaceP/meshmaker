""" Shapefile to .Geo file converter

Module to convert a shapefile '.shp' to a geometry file of gmsh '.geo'.
!!! IMPORTANT !!! It can handle inner polygons BUT it can't place boundaries in the inner polygons!

author: Ignace Pelckmans
            (University of Antwerp, Belgium)
"""

import pickle
import time
import warnings
warnings.filterwarnings("ignore")
import numpy as np
import matplotlib.pyplot as plt

from shp2mpol import *
from pol2pol import *
from remove_vertices import *
from remove_duplicates import *
from polplot import *
from np2bgfield import *

#############################################################################################################
# Parameter Definition
#############################################################################################################

# --- Shapefile Parameters --- #
remove_vertices = True # remove vertices if there is no change of direction
remove_duplicates = True # remove duplicate coordinate pairs, highly recommended!
boundaries = [[1,2]]


# --- Plotting Flags --- #
plot_flag = True
plot_vert = True
plot_vertlabels = True
plot_each_x_labels = 5
figfile = '/home/ignace/Desktop/test.jpg'
figure_size = (50,50)

# --- Take interior rings of shapefile into account? --- #
inner_flag = True
include_inner = False

# --- Generate a .geo file? --- #
generate_geo_file = True

# --- Spline parameters --- #
spline_max_len = 50 # max number of vertices in a spline, if no splits are desired: spline_max_len =  float('inf')
breaks = [] # vertices at which to split the splines

# --- Input and output directories --- #
inputfile = '/home/ignace/Documents/PhD/Data/Land Use/Guayas_Land_Use/Churute_Subregion_ROI.shp'
proj = 32717 # epsg code of the desired projection (hint: use metric projections, for instane, avoid epsg:4326)
outputfile_geo = '/home/ignace/Documents/PhD/TELEMAC/Playing_Around/First_Steps/Test_Cases/input/v1/Churute_Subregion.geo'
outputfile_mesh = '/home/ignace/Documents/PhD/TELEMAC/Playing_Around/First_Steps/Test_Cases/input/v1/Churute_Subregion.msh'



# --- Mesh Properties --- #
heterogenuous = True
homo_mesh_size = 250 # if heterogenuous is False and thus a homogeneous mesh is requested
XY_field = './Churute_subregion/bf_params.pck' # if heterogenuous is True
# path to a pickle which hold a list with the following elements:
#    - 2d numpy array representing a spatial raster with values for the mesh size
#    - coordinate pair with the coordinates of the left bottom corner of the 2d numpy array
#    - x resolution
#    - y resolution

outside_mesh_size = 200 # mesh size outside the given mesh size field
multiplier = 0.33
minimum_resolution = 5
interpolate_buffer_thickness = 350

#############################################################################################################

# track time
t0 = time.time()

#-----------------------------------------------#
#-- Load Polygon and prepare coordinate pairs --#
#-----------------------------------------------#

# load .shp and reproject
pol, epsg = shp2Mpol(inputfile, return_coordinate_system = True)
pol = pol2Pol(pol, epsg, proj)[0]

# remove excessive vertices
if remove_vertices: pol = removeVertices(pol)[0]
# Remove duplicates (highly recommended!)
if remove_duplicates: pol = removeDuplicates(pol)[0]

# convert to Numpy Array and rotate 90° do the shape becoms n x 2, remove z-coordinates
exter = np.rot90(pol.exterior.xy)[:-1,0:2]
inter = [np.rot90(i.xy)[:-1,0:2] for i in pol.interiors]

# number of exterior vertices
n = len(exter)

#-----------------------------------------------#
#--------- Plot Polygon, if indicated ----------#
#-----------------------------------------------#

# if indicated, plot the original polygon and possibly (communicated through flags) vertices and their labels.
# Do this to find the id's of the boundaries and if needed, where to place breaks

if plot_flag:
    f, a = plt.subplots(figsize = figure_size)
    a.set_aspect('equal')
    # plot polygon
    if inner_flag: polPlot(exter, XY_inner = inter, show_vertices = plot_vert, empty= False, plot_on_axis = a,
                               show_vertices_labels= plot_vertlabels, show_vertices_labels_interval = plot_each_x_labels, vertices_color = 'maroon')
    else: polPlot(exter, show_vertices = plot_vert, empty= False, plot_on_axis = a,
                      show_vertices_labels= plot_vertlabels, show_vertices_labels_interval = plot_each_x_labels, vertices_color = 'maroon')

    # indicate boundaries
    for b in boundaries:
        X,Y = [[],[]]
        for i in range(b[0], b[1]+1):
            x, y = exter[i-1]
            X.append(x)
            Y.append(y)
        bound = a.plot(X, Y, color='purple', label = 'boundaries')
    # indicate break points
    for b in breaks:
        x, y = exter[b-1]
        bre = a.scatter(x,y, c = 'pink', label = 'breaks')

    f.savefig(figfile)


#-----------------------------------------------#
#------------ Write the .geo file --------------#
#-----------------------------------------------#

# only proceed if we want to generate a geo file
if generate_geo_file:

    # Open txt file and allow both writing and reading
    f = open(outputfile_geo, "w+")

    # write header to indicate the right engine
    f.write("SetFactory(\"OpenCASCADE\");\n//+\n")

    # initialize counters for each features
    n_Point = 1
    n_Spline = 1
    n_LineLoop = 1
    n_Surface = 1

    # reorganize XY list to put boundary first
    # -----------------------------------------------------------------

    # there is a spline-break at the first vertex anyway, so if there is a boundary, place the first vertex of that boundary
    # as the first vertex so the break will be at a boundary. To do so, we go over all vertices id's and store them in a list
    # when we pass the vertex id which is the first vertex of the first boundary, we start storing the vertex id's in a new list.
    # The first list will now be stiched after the second one.
    #
    # For example, if we have an exterior of 10 vertices and there is one boundary between 4 - 5:
    # original exterior vertex id's list: [1, 2, 3, 4, 5, 6, 7, 8, 9, 10]
    # new exterior vertex id's list: [4, 5, 6, 7, 8, 10, 1, 2, 3]

    # only proceed if there are boundaries
    if len(boundaries)>0:

        # id of the first vertex of the first boundary
        b0 = boundaries[0][0]

        # initiate two lists
        exter_re1, exter_re2 = [], []

        # loop over all vertices in the exterior
        for i in range(len(exter)):
            # if i + i (vertex id's start at 1) is smaller than the first vertex of the first boundary ...
            if i + 1 < b0:
                # ... append that vertex to the first list
                exter_re1.append(exter[i])
            else:
                # otherwise, append to the other list
                exter_re2.append(exter[i])

        # rearange order by placing the second list first
        exter_re2.extend(exter_re1)
        # exter is now the updated exter list
        exter = exter_re2

        # this means that the id's are reorganized and so, also the boundaries and breaks should be reorganized
        boundaries_upd = []
        # reorganize all boundary id's
        #
        # example: [1, 2, 3, 4, 5, 6, 7, 8, 9, 10] with boundaries at 4 - 5 and 8 - 10
        # Here, the new boundaries will be: 1 - 2 and 5 - 7

        for b_l, b_r in boundaries:
            boundaries_upd.append([b_l - b0 + 1, b_r - b0+1])
        boundaries = boundaries_upd

        # similar for all break points
        # however, now break id's can be smaller than 0 and thus - b0 + 1 can result in updated breaks equal to 0 or 1
        # in that case, add the total number of exterior breaks
        #
        # example: [1, 2, 3, 4, 5, 6, 7, 8, 9, 10]
        #              *     ^''^  *     ^'''''^
        # with boundaries at 4 - 5 and 8 - 10 and breaks at 2 and 6
        #
        # updated example: [1, 2, 3, 4, 5, 6, 7, 8, 9, 10]
        #                   ^''^  *     ^'''''^     *
        # Here, the new boundaries will be: 1 - 2 and 5 - 7
        # and the new breaks will be: 3
        breaks_upd = []

        breaks = [b-b0+1 for b in breaks]
        breaks = [b+n if b < 1 else b for b in breaks ]

    # Add Points
    # -----------------------------------------------------------------

    # add a header
    f.write('//+\n//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++\n'
            '//+ Points\n//++++++++++++++++++++++++++++++++++++++++++++++++++++'
            '++++++++++++++++++\n')

    # each exterior vertex should be added as a Point in the format:
    # Point(1) = {x, y, 0}
    # Point(2) = {x, y, 0}
    # ...

    for i in range(np.shape(exter)[0]):
            f.write("Point(%d) = {%f, %f, 0};\n" % (n_Point, exter[i][0], exter[i][1]))
            # keep track of the number of points
            n_Point += 1


    # Add Splines and Line
    # -----------------------------------------------------------------

    # add a header
    f.write('//+\n//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++\n'
            '//+ Lines & Splines\n//++++++++++++++++++++++++++++++++++++++++++++++++++++'
            '++++++++++++++++++\n')

    # if there are boundaries
    if len(boundaries) > 0:
        # initialize list to store the id's of the boundaries (or 'physical lines')
        physical_lines = []

        # start with first boundary since the first vertices will always be boundary vertices in case of boundaries

        # first boundary
        bound1 = boundaries[0]
        # length of first boundary
        bound1diff = bound1[1]-bound1[0]

        # write a line object in the .geofile
        # format example:
        # Line(1) = {1, 2}
        # where the numbers inbetween curly brackets refer to .geo point id's
        f.write('Line(%d) = {1, ' % n_Spline)

        # add this id to the list of physical boundaries
        physical_lines.append(n_Spline)

        # add a line per segment of the first boundary
        # for instance if there is a boundary 1 - 3:
        # Line(1) = {1, 2}
        # Line(2) = {2, 3}
        for i in range(0,bound1diff-1):
            f.write('%d};\n' % (i+2))
            # keep track of number of lines and splines
            n_Spline += 1
            f.write('Line(%d) = {%d,' % (n_Spline, i+2))
            # add all these lines to the list of physical boundaries
            physical_lines.append(n_Spline)

        # finish the open Line definition
        f.write('%d};\n' % (bound1diff+1))
        n_Spline += 1
        # write the last line of the first boundary
        f.write('Spline(%d) = {%d, ' % (n_Spline, bound1diff+1))
        n_Spline += 1

        # variable to keep track if the loop is inside a boundary or not
        b = False
        # initiate a variable to count the number of segments in a spline
        spline_counter = 0
        # loop over all points except for the points in the first boundary
        for i in range(bound1diff+2,n_Point):
            # if the point is the starting vertex of a boundary, close whatever was open en start a new line
            # Example:
            #
            # before: Spline(n) = {1, 2, 3 .
            # newly added:
            # , 4};
            # Line(n + 1) = Line{4,
            #
            if list(np.asarray(boundaries)[:,0]).__contains__(i):
                # start a new line
                f.write("%d};\nLine(%d) = {%d, " % (i, n_Spline, i))
                physical_lines.append(n_Spline)
                n_Spline += 1
                # the loop is inside a boundary
                b = True

                # if the boundary only exist out of one vertex
                if list(np.asarray(boundaries)[:, 1]).__contains__(i+1):
                    b = False

            # if the vertex is part of a boundary
            elif b:
                # close the upper line and start a new one
                # Example:
                #
                # before:
                # Line(n + 1) = Line{4,
                # new:
                # Line(n + 1) = Line{4, 5};
                # Line(n + 2) = Line{5,
                f.write('%d};\nLine(%d) = {%d,' % (i, n_Spline, i))
                physical_lines.append(n_Spline)
                n_Spline += 1
                # close the boundary flag if this vertex is the second last vertex of a boundary
                if list(np.asarray(boundaries)[:, 1]).__contains__(i+1):
                    b = False

            # if the vertex is the closing vertex of a boundary
            elif list(np.asarray(boundaries)[:, 1]).__contains__(i):
                # close the upper line and start a new one
                # Example:
                #
                # before:
                # Line(n + 2) = Line{5,
                # new:
                # Line(n + 2) = Line{5, 6};
                # Line(n + 3) = Line{6,
                f.write("%d};\nSpline(%d) = {%d," % (i, n_Spline,i))
                n_Spline += 1
                # if this point is also the last point
                # Example:
                #
                # before:
                # Line(n + 3) = Line{6,
                # new:
                # Line(n + 3) = Line{6, 7};

                if i == n_Point-1:
                    f.write("%d, 1};\n" % (n_Point-1))


            # if this point is also the last point (and no part of a boundary)
            # Example:
            #
            # before:
            # Line(n + 3) = Line{6,
            # new:
            # Line(n + 3) = Line{6, 7};
            elif i == n_Point-1:
                f.write("%d, 1};\n" % (n_Point-1))

            # if the vertex is just a regular vertex and no part of a boundary
            else:
                # count the number of vertices in a spline
                spline_counter += 1
                # if there is a max spline length indicated
                if spline_max_len:
                    # if the number of vertices in the spline is lower than than the allowed number
                    if spline_counter <= spline_max_len:
                        # if the vertex is a break
                        if breaks.__contains__(i):
                            # close the upper spline and start a new one
                            f.write("%d};\nSpline(%d) = {%d, " % (i, n_Spline, i))
                            n_Spline += 1
                        else:
                            # if the vertex is not a break and the number of vertices is lower than the allowed number, just add the point to the open spline
                            f.write("%d, " % (i))
                    else:
                        # if the number of vertices in this spline would exceed the allowed number, start a new spline
                        f.write("%d};\nSpline(%d) = {%d, " % (i, n_Spline, i))
                        n_Spline += 1
                        spline_counter = 0
                # if there is not a maximum number of vertices in a spline indicated, check if it is a break vertex
                elif breaks.__contains__(i):
                        f.write("%d};\nSpline(%d) = {%d, " % (i, n_Spline, i))
                        n_Spline += 1
                # if there is not a maximum number of vertices in a spline and if it is not a break vertex, just add it to the open spline
                else:
                    f.write("%d, " % (i))

    # if there are no boundaries indicates
    else:
        # create first spline
        f.write('Spline(1) = {')
        # add all points to this spline, except for breaks. In that case, start a new spline
        for i in range(1,n_Point):
            # in case of breaker vertices
            if breaks.__contains__(i):
                    f.write("%d};\nSpline(%d) = {%d, " % (i, n_Spline, i))
                    n_Spline += 1
            f.write('%d,' % (i))
        f.write('%d};\n' % (n_Point))


    # Inner holes
    #------------------------------------------------------------------
    # if there are inner rings which should be 'islands' in the domain
    if inner_flag:
        # check if there are inner rings
        if len(inter)>0:
            # write a header
            f.write('//+\n//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++\n'
                    '//+ Holes Points\n//++++++++++++++++++++++++++++++++++++++++++++++++++++'
                    '++++++++++++++++++\n')

            # save all inner points in a list per ring
            Inner_Points = []
            # loop over all inter rings
            for I in inter:
                Inner_Points.append([])
                # loop over all vertices in an inner ring
                for i in range(0,len(I)):
                    # add as a point
                    f.write("Point(%d) = {%f, %f, 0};\n" % (n_Point, I[i][0], I[i][1]))
                    # save the point id's in a list of list per inner ring
                    Inner_Points[-1].append(n_Point)
                    n_Point += 1

            # start a new counter to include all inner rings
            n_Spline_inter = n_Spline
        if len(inter)>0:
            # create header
            f.write('//+\n//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++\n'
                    '//+ Holes Splines\n//++++++++++++++++++++++++++++++++++++++++++++++++++++'
                    '++++++++++++++++++\n')

            # initiate a list to store all inner spline ids
            Inner_Splines_ids = []
            # create a spline per inner ring (no breaks or boundaries are allowed at this point)
            for I in Inner_Points:
                # open a spline
                f.write('Spline(%d) = {' % (n_Spline_inter))
                # add the id to the list of inner splines
                Inner_Splines_ids.append(n_Spline_inter)
                # count
                n_Spline_inter += 1
                # add all points to this spline
                for i in I:
                    f.write('%d,' % (i))
                # end with the first point of the ring
                f.write('%d};\n' % (I[0]))


    # add physical boundaries
    #------------------------------------------------------------------3
    # add a header
    f.write('//+\n//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++\n'
            '//+ Boundaries\n//++++++++++++++++++++++++++++++++++++++++++++++++++++'
            '++++++++++++++++++\n')
    # if there are no boundaries
    if len(boundaries)==0:
        for i in range(n_Spline):
            f.write('Physical Line("closed") = {%d};\n' % (i+1))
    # if there are boundaries
    else:
        # initialize a list to store all id's of what should be closed lines
        closed_lines = []
        # loop over all
        for i in range(1, n_Spline):
            # if it is not a physical line
            if not physical_lines.__contains__(i):
                closed_lines.append(i)
    # add all boundary lines to the physical line list with label "Tidal Entrance"
    if len(boundaries) > 0:
        f.write('Physical Line("TidalEntrance") = {')
        for i in range(len(physical_lines)-1):
            # add each pysical line
            f.write('%d,' % (physical_lines[i]))
        # close the physical line object
        f.write('%d};\n' % (physical_lines[-1]))

    # non-boundary lines are 'closed' boundaries
    f.write('Physical Line("closed") = {')
    for i in range(len(closed_lines)-1):
        f.write('%d,' % (closed_lines[i]))
    f.write('%d};\n' % (closed_lines[-1]))

    # all inner holes are closedInner boundaries
    f.write('Physical Line("closedInner") = {')
    for i in range(len(Inner_Splines_ids)-1):
        f.write('%d,' % (Inner_Splines_ids[i]))
    f.write('%d};\n' % (Inner_Splines_ids[-1]))

    # add two lines to define the a priori resolution
    #------------------------------------------------------------------

    # only needed if the mesh size is set at homogeneous
    if not heterogenuous:
        # add header
        f.write('//+\n//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++\n'
                '//+ A Priori Resolution\n//++++++++++++++++++++++++++++++++++++++++++++++++++++'
                '++++++++++++++++++\n')
        # Accuracy of evaluation of the LC field for 1D mesh generation
        # Prohibit boundary mesh size interpolation (suggested if the cell size is indicated by a field)
        f.write('Mesh.LcIntegrationPrecision = 1e-3;\nMesh.CharacteristicLengthExtendFromBoundary = 0;\n')

    # Create line loop and Surface
    # -----------------------------------------------------------------

    # create a header
    f.write('//+\n//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++\n'
            '//+ Create Line Loop and Surface\n//++++++++++++++++++++++++++++++++++++++++++++++++++++'
            '++++++++++++++++++\n')

    # create a line loop with all splines of the exterior
    # open the line loop
    f.write('Line Loop(%d) = {' % (n_LineLoop))
    # add all the splines
    for i in range(1,n_Spline-1):
        f.write('%d, '% (i))
    # close the line loop
    f.write('%d};\n' % (n_Spline-1))
    # create a plane surface from the lineloop
    f.write('Plane Surface(%d) = {%d};\n' % (n_Surface, n_LineLoop))
    # create a physical surface from the plane surface (no water can pass this surface)
    f.write('Physical Surface("inside") = {%d};\n' % (n_Surface))

    # count!
    n_LineLoop += 1
    n_Surface += 1

    # Inner holes - Line loop and Surfaces
    #------------------------------------------------------------------

    # check if there are inner rings
    if len(inter) > 0:
        # create header
        f.write('//+\n//+ Inner loops and surfaces\n//+\n')
        # counter
        t = 0
        # go over all inner rings ( Inner_Points is a list of list and so, its length is equal to the number of rings)
        for I in Inner_Points:
            if include_inner == False:
                # per ring, write a line loop, plane and physical surface
                f.write('Line Loop(%d) = {%d};\n' % (n_LineLoop, n_Spline + t))
                f.write('Plane Surface(%d) = {%d};\n' % (n_Surface, n_LineLoop))
                # remove this inner ring from the outer ring
                f.write('BooleanDifference(%d) = { Surface{%d}; Delete; }{ Surface{%d}; Delete;};\n' % (n_Surface+1, n_Surface-1, n_Surface))
            elif include_inner == True:
                # per ring force that spline on the surface
                f.write('Line{%d} In Surface{%d};\n' % (n_Spline + t, n_Surface -t*2))
            # count!
            t += 1
            n_LineLoop += 1
            n_Surface += 2


        # create a final physical surface to set the domain
        f.write('Physical Surface("domain") = {%d};\n' % (n_Surface - 1))

    # Set Homogenuous mesh
    # -----------------------------------------------------------------
    # if  homogeneous
    if not heterogenuous:
        # add header
        f.write('//+\n//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++\n'
                '//+ Set Homogeneous Mesh\n//++++++++++++++++++++++++++++++++++++++++++++++++++++'
                '++++++++++++++++++\n')
        # create a field with the set homogeneous cell size everywhere
        f.write('//+\nField[1] = MathEval;\nField[1].F = "%f";\n' % (homo_mesh_size))
        # set this field as the background field
        f.write('Background Field = 1;')
    # if the cell size is set as heterogenuous
    else:
        # open the pickled file with the needed information for a heterogenuous background field
        with open(XY_field, 'rb') as input:
            xy = pickle.load(input)
        # background file
        background_file = outputfile_geo[:-3] + 'bg'
        # function to create a valid background field from a Numpy 2D array
        np2BackgroundField(xy[0], background_file , multiplier = multiplier, minres = minimum_resolution,
            zerovalues = outside_mesh_size, zeropoint= xy[1], xres = xy[2], yres = xy[3], plot_arr = './test.png',
            buffer_interpolation = interpolate_buffer_thickness)

        # add header
        f.write('//+\n//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++\n'
                '//+ Set Heterogenuous Mesh\n//++++++++++++++++++++++++++++++++++++++++++++++++++++'
                '++++++++++++++++++\n')

        # complete the background directory if it's a local directory
        if background_file[0] == '.':
            background_file = os.getcwd() + background_file[1:]
        # add this background field as the structured background field
        f.write('Field[1] = Structured;\nField[1].FileName = "%s";\n' % (background_file))
        # set the cell size outside the background field
        f.write('Field[1].OutsideValue = %f;\n' % (outside_mesh_size))
        f.write('Field[1].SetOutsideValue = 1;\n')
        # set the format of this background field and set is at the background file to apply the meshing on
        f.write('Field[1].TextFormat = 1;\nBackground Field = 1;\n//+')

    # close the .geo file
    f.close()

    # mesh the .geo file
    os.system('gmsh ' + outputfile_geo + ' -2 -o ' + outputfile_mesh + ' -format msh2')

    # show how long it took to generate this file
    t = time.time() - t0
    min = t//60
    sec = t%60
    print('Ready! .geo-file is written in %s. It took %.2f minutes and %.2f seconds.' % (outputfile_mesh, min, sec))
