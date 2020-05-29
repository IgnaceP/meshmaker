""" remove Vertices

Function to remove excessive vertices in a polygon
Excessive is defined as vertex which are not located on a corner and thus, the direction
of the segment does not change.

author: Ignace Pelckmans
    (University of Antwerp, Belgium)
"""

import numpy as np
from shapely.geometry import MultiPolygon, Polygon

def removeVerticesFromArray(xy, geotype = 'pol'):
    """
     Function to remove excessive vertices of line/polygon
     Excessive is defined as vertex which are not located on a corner and thus, the direction
     of the segment does not change.

     Args:
        xy: (Required) Numpy array with dimensions (n x 2) or a list of pairs representing x- and y-coordinates of the polygon's vertices
        geotype: (Optional, defaults to 'pol')'pol' or 'line' indicating if it is a line or a closed polygon

    Returns:
      Numpy array with dimensions (n x 2)  representing x- and y-coordinates of the polygon's vertices without the excessive points
    """

    # Make sure it can handle a list or a Numpy array
    if type(xy) == list:
        xy = np.asarray(xy)

    # get number of vertices
    n = np.shape(xy)[0]

    # initialize an empty numpy array to store the vertices to be removed
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

        # check slope of straights formed between i and it left/right neigbor
        x, y = xy[i,:]
        lx, ly = xy[l,:]
        rx, ry = xy[r,:]
        la = (ly - y)/(lx - x)
        ra = (ry - y)/(rx - x)

        # if slopes are equal, the vertex is not essential
        if abs(la - ra) <= 0:
            rem[i] = 1

    # mask to remove all excessive vertices
    xy_upd = xy[rem==0]

    return xy_upd

def removeVertices(mpol, geotype = 'pol'):
    """
     Function to remove excessive vertices of shapely (Multi)Polygon
     Excessive is defined as vertex which are not located on a corner and thus, the direction
     of the segment does not change.

     Args:
        xy: (Required) shapely (Multi)Polygon to edit
        geotype: (Optional, defaults to 'pol')'pol' or 'line' indicating if it is a line or a closed polygon

    Returns:
      Shapely (Multi)Polgyon without the excessive vertices
    """

    # check if it is a Polygon or a MultiPolygon
    if type(mpol) == Polygon:
        mpol = [mpol]

    # Initialize list for new pols
    pols = []

    # Loop over all Polygons
    for pol in mpol:

        # exterior
        xy = np.rot90(pol.exterior.xy)
        xy = removeVerticesFromArray(xy)

        # interiors
        xys_i = []
        for i in pol.interiors:
            xy_i = np.rot90(i.xy)
            xy_i = removeVerticesFromArray(xy_i)
            xys_i.append(xy_i)

        # create Polygon and add to polygon list
        pols.append(Polygon(xy, xys_i))

    # if the input was a Polygon, return a Polgyon
    if type(mpol) == Polygon:
        return pols[0]
    # if the input was a MultiPolygon, return a MultiPolygon
    else:
        return MultiPolygon(pols)
