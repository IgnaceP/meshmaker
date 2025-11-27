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
        mpol = MultiPolygon([mpol])

    mpol = mpol.simplify(0)

    # if the input was a Polygon, return a Polgyon
    if type(mpol) == Polygon:
        return pols[0]
    # if the input was a MultiPolygon, return a MultiPolygon
    else:
        return MultiPolygon(pols)
