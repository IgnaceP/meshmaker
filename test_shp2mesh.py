from shp2mesh import pymesh

HOME = '/home/ignace/Documents/PhD/TELEMAC/Playing_Around/First_Steps/Test_Cases/input/v3/'

msh = pymesh('test')
msh.loadShp('/home/ignace/Documents/PhD/Data/Land Use/Guayas_Land_Use/Churute_Subregion_ROI.shp')
msh.setBoundariesAndBreaks(boundaries =[[1,2]])
msh.plotShp('./testPlotpol.png', plot_inner = False, figure_size = (10,10), plot_vertlabels = True, plot_each_x_labels = 1)
msh.writeGeoFile(HOME + 'Churute_Subregion.geo',
                cell_size = 'hetero',
                inner_rings = True,
                XY_field = HOME + 'bf_params.pck',
                multiplier = 0.33,
                minres = 25,
                zerovalues = 200,
                outside_mesh_size = 200,
                plot_arr = HOME + 'Churute_Subregion_backgroundfield.png',
                buffer_interpolation = 200,
                algo = 'Delaunay',
                physical_line_labels = ['closed', 'TidalEntrance', 'OpenInner'])

msh.meshTheDomain(HOME + 'Churute_Subregion.msh')
inner_holes_xy = msh.getCoordinatesInnerHoleNodes(HOME + 'Churute_Subregion.msh')
msh.meshInnerHoles(HOME + 'Churute_Subregion_incl.msh', inner_holes_xy)
