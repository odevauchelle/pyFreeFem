import matplotlib.pyplot as pp
from tempfile import NamedTemporaryFile
from random import sample

import sys
sys.path.append('./../')

import pyFreeFem as pyff

# Create mesh with FreeFem++

edp_str = '''
border Circle( t = 0, 2*pi ){ x = cos(t); y = sin(t); }
mesh Th = buildmesh( Circle(150) );
'''

edp_str += pyff.export_mesh_edp() # adds a few lines to edp string

FreeFem_output = pyff.run_FreeFem( edp_str )

mesh = pyff.FreeFem_str_to_mesh( FreeFem_output )

# Change mesh

for triangle_index in sample( range( len(mesh.triangles ) ), 20 ) :
    mesh.boundary_edges.update( { ( triangle_index, 0 ) : 2 } )

# Rename boundaries
new_names = {1:'initial', 2:'new'}
mesh.rename_boundary( new_names )


# Export mesh back to FreeFem

temp_mesh_file = NamedTemporaryFile( suffix = '.msh' )
pyff.savemesh( filename = temp_mesh_file.name, mesh = mesh )

edp_str = '''
mesh Th = readmesh( "mesh_file_name" ) ;
'''.replace( 'mesh_file_name', temp_mesh_file.name )

edp_str += pyff.export_mesh_edp()

# calculate FEM matrices

edp_str +='''
fespace Vh( Th, P1 ) ;
Vh u, v ;
'''

matrix_types = [ pyff.stiffness, pyff.Grammian, pyff.boundary_Grammian( 1, 2 ) ]

for matrix_type in matrix_types :
    edp_str += pyff.export_matrix_edp( **matrix_type )

FreeFem_output = pyff.run_FreeFem( edp_str )

mesh = pyff.FreeFem_str_to_mesh( FreeFem_output )
mesh.rename_boundary( {1:'initial', 2:'new'} )

matrices = {}

for matrix_type in matrix_types :
    matrices[ matrix_type['matrix_name'] ] = pyff.FreeFem_str_to_matrix( FreeFem_output, matrix_type['matrix_name'] )


##########################################################################################################

from scipy.sparse.linalg import spsolve
import numpy as np

epsilon = 1e-4
M = - matrices[ pyff.stiffness['matrix_name'] ] + 1./epsilon*matrices[ pyff.boundary_Grammian(1,2)['matrix_name'] ]
Source = matrices[ pyff.Grammian['matrix_name'] ]*np.array( [1]*len( mesh.x ) )
u = spsolve( M, Source )


##################
#
# First FIGURE
#
#################

figs = {}

figs['mesh'] = pp.figure()
mesh.plot_triangles( color = 'k', lw = .5, alpha = .2  )
mesh.plot_boundaries()
pp.legend(title = 'Boundary')
x_lim, y_lim = pp.gca().get_xlim(), pp.gca().get_ylim()

figs['field'] = pp.figure()
pp.tricontourf( mesh, u )
mesh.plot_boundaries( color = 'black' )


for fig_key in figs.keys() :

    pp.figure(figs[fig_key].number)

    pp.axis('equal')
    pp.axis('off')
    pp.xticks([])
    pp.yticks([])

    pp.xlim(x_lim); pp.ylim(y_lim)

    pp.savefig( '../figures/' + __file__.split('/')[-1].split('.')[0] + '_' + fig_key + '.svg' , bbox_inches = 'tight' )





pp.show()
