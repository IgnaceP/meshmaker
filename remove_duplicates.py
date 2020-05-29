""" remove Duplicates

Remove all duplicate vertices from a shapely (MultiPolygon)

author: Ignace Pelckmans
    (University of Antwerp, Belgium)
"""

import numpy as np
from shapely.geometry import MultiPolygon, Polygon

def removeDuplicates(mpol):
    """
     Function to remove duplicate vertices from a shapely (Multi)Polygon

     Args:
        mpol: (Required) shapely (Multi)Polygon

    Returns:
      Numpy array with dimensions (n x 2)  representing x- and y-coordinates of the polygon's vertices without the excessive points
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
        indexes = np.unique(xy, return_index=True, axis = 0)[1]
        xy = np.asarray([xy[index] for index in sorted(indexes)])

        # interiors
        xys_i = []
        for i in pol.interiors:
            xy_i = np.rot90(i.xy)
            indexes = np.unique(xy_i, return_index=True, axis = 0)[1]
            xy_i = np.asarray([xy_i[index] for index in sorted(indexes)])
            xys_i.append(xy_i)

        # create Polygon and add to polygon list
        pols.append(Polygon(xy, xys_i))

    # if the input was a Polygon, return a Polgyon
    if type(mpol) == Polygon:
        return pols[0]
    # if the input was a MultiPolygon, return a MultiPolygon
    else:
        return MultiPolygon(pols)
