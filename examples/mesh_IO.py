import sys
sys.path.append('./../')

import pyFreeFem as pyff
import matplotlib.pyplot as pp
from tempfile import NamedTemporaryFile


##################
#
# CREATE MESH
#
#################

edp_str = '''
border Circle( t = 0, 2*pi ){ x = cos(t); y = sin(t); }
mesh Th = buildmesh( Circle(10) );
'''

edp_str += pyff.export_mesh_edp() # adds a few lines to edp string

FreeFem_output = pyff.run_FreeFem( edp_str )

mesh = pyff.FreeFem_str_to_mesh( FreeFem_output )

segments = []

for triangle in mesh.triangles[2:3]:
    print(triangle[1:])

    segments += [ [triangle[1:]] ]

mesh.boundaries += [ pyff.Boundary( segments, label = '12' ) ]

##################
#
# FIRST FIGURE
#
#################
pp.figure()
mesh.plot_triangles(indices = True)
mesh.plot_nodes(indices = True)
mesh.plot_boundaries()

pp.axis('equal')
pp.axis('off')
pp.xticks([])
pp.yticks([])

pp.savefig( '../figures/' + __file__.split('/')[-1].split('.')[0] + '_1.svg' , bbox_inches = 'tight' )

##################
#
# EXPORT MESH TO FREEFEM
#
#################


temp_mesh_file = NamedTemporaryFile( suffix = '.msh' )
mesh.save( filename = temp_mesh_file.name )

edp_str = '''
mesh Th = readmesh( "mesh_file_name" ) ;
'''.replace( 'mesh_file_name', temp_mesh_file.name )

edp_str += pyff.export_mesh_edp() # adds a few lines to edp string

FreeFem_output = pyff.run_FreeFem( edp_str )

mesh = pyff.FreeFem_str_to_mesh( FreeFem_output )

##################
#
# SECOND FIGURE
#
#################

pp.figure()
mesh.plot_triangles(indices = True)
mesh.plot_nodes(indices = True)

mesh.plot_boundaries()

pp.axis('equal')
pp.axis('off')
pp.xticks([])
pp.yticks([])

pp.savefig( '../figures/' + __file__.split('/')[-1].split('.')[0] + '_2.svg' , bbox_inches = 'tight' )

pp.show()
