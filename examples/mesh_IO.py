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
mesh Th = buildmesh( Circle(20) );
'''

edp_str += pyff.export_mesh_edp() # adds a few lines to edp string

FreeFem_output = pyff.run_FreeFem( edp_str )

mesh = pyff.FreeFem_str_to_mesh( FreeFem_output )

mesh.boundaries += [ pyff.Boundary( [[ mesh.triangles[3][1:] ]], label = 'B' ) ]

##################
#
# FIRST FIGURE
#
#################

mesh.plot_triangles( labels = True )
mesh.plot_nodes()
mesh.plot_boundaries( labels = True )

pp.axis('equal')
pp.axis('off')
pp.xticks([])
pp.yticks([])

pp.savefig( '../figures/' + __file__.split('/')[-1].split('.')[0] + '_1.svg' , bbox_inches = 'tight' )

pp.show()

##################
#
# EXPORT MESH TO FREEFEM
#
#################


temp_mesh_file = NamedTemporaryFile( suffix = '.msh' )

edp_str = '''
mesh Th = readmesh( "mesh_file_name" ) ;
'''.replace( 'mesh_file_name', temp_mesh_file.name )

edp_str += pyff.export_mesh_edp() # adds a few lines to edp string

print(edp_str)

FreeFem_output = pyff.run_FreeFem( edp_str )

mesh = pyff.FreeFem_str_to_mesh( FreeFem_output )

##################
#
# SECOND FIGURE
#
#################

mesh.plot_triangles( labels = True )
mesh.plot_nodes()
mesh.plot_boundaries( labels = True )

pp.axis('equal')
pp.axis('off')
pp.xticks([])
pp.yticks([])

pp.savefig( '../figures/' + __file__.split('/')[-1].split('.')[0] + '_2.svg' , bbox_inches = 'tight' )

pp.show()
