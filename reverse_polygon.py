from shapely.geometry import Polygon

def reversePolygon(geometry):
    """Return a new geometry with the reversed order of vertices in all of its linear rings."""

    # Define a recursive function to reverse the order of vertices in a linear ring
    def reverse_linear_ring(ring):
        return type(ring)(list(ring.coords)[::-1])

    if isinstance(geometry, Polygon):
        # Reverse the order of vertices in the exterior and interior rings of the polygon
        reversed_exterior = reverse_linear_ring(geometry.exterior)
        reversed_interiors = [reverse_linear_ring(interior) for interior in geometry.interiors]

        # Create a new polygon with the reversed rings
        reversed_geometry = Polygon(shell=reversed_exterior, holes=reversed_interiors)
    elif isinstance(geometry, MultiPolygon):
        # Reverse the order of vertices in each polygon of the multipolygon
        reversed_geometries = [reverse_geometry(polygon) for polygon in geometry]

        # Create a new multipolygon with the reversed polygons
        reversed_geometry = MultiPolygon(reversed_geometries)
    else:
        raise ValueError("Input geometry must be a Polygon or MultiPolygon object.")

    return reversed_geometry

