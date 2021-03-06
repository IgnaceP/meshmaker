from shp2mesh import pymesh

HOME = '/home/ignace/Documents/PhD/TELEMAC/Playing_Around/First_Steps/Test_Cases/input/v2/'

msh = pymesh('test')
msh.loadShp('/home/ignace/Documents/PhD/Data/Land Use/Guayas_Land_Use/Churute_Subregion_ROI.shp')
msh.setBoundariesAndBreaks(boundaries =[[4,6], [11,13]])
msh.plotShp('./test.png', plot_inner = False, figure_size = (10,10))
msh.writeGeoFile(HOME + 'Churute_Subregion.geo',
                cell_size = 'hetero',
                inner_rings = True,
                include_inner_rings = False,
                XY_field = HOME + 'bf_params.pck',
                multiplier = 0.33,
                minres = 25,
                zerovalues = 200,
                outside_mesh_size = 200,
                plot_arr = HOME + 'Churute_Subregion_backgroundfield.png',
                buffer_interpolation = 200)

msh.meshTheDomain('/home/ignace/Documents/PhD/TELEMAC/Playing_Around/First_Steps/Test_Cases/input/v2/Churute_Subregion.msh')
