import sys
sys.path.append('./../')

import pyFreeFem as pyff
import matplotlib.pyplot as pp

edp_str = '''
real smallRadius = .3;
border outerCircle( t = 0, 2*pi ){ x = cos(t); y = 0.8*sin(t); }
border innerCircle( t = 2*pi, 0 ){ x = .5 + smallRadius*cos(t); y = smallRadius*sin(t); }
mesh Th = buildmesh( outerCircle(100) + innerCircle(40) );

fespace Vh( Th, P1 );
Vh u,v;
'''

edp_str += pyff.export_mesh_edp()

matrix_types = [ pyff.stiffness, pyff.Grammian, pyff.boundary_Grammian(1,2) ]

for matrix_type in matrix_types :
    edp_str += pyff.export_matrix_edp( **matrix_type )

FreeFem_output = pyff.run_FreeFem( edp_str )

mesh = pyff.FreeFem_str_to_mesh( FreeFem_output, simple_boundaries = True )

matrices = {}

for matrix_type in matrix_types :
    matrices[ matrix_type['matrix_name'] ] = pyff.FreeFem_str_to_matrix( FreeFem_output, matrix_type['matrix_name'] )

mesh.plot_triangles( lw = 1 )
mesh.plot_boundaries( labels = False )

##########################################################################################################

pp.axis('equal')
pp.axis('off')
pp.xticks([])
pp.yticks([])

pp.savefig( '../figures/' + __file__.split('/')[-1].split('.')[0] + '.svg' , bbox_inches = 'tight' )

x_lim, y_lim = pp.gca().get_xlim(), pp.gca().get_ylim()

pp.show()

##########################################################################################################

from scipy.sparse.linalg import spsolve
import numpy as np

epsilon = 1e-4
M = - matrices[ pyff.stiffness['matrix_name'] ] + 1./epsilon*matrices[ pyff.boundary_Grammian(1,2)['matrix_name'] ]
Source = matrices[ pyff.Grammian['matrix_name'] ]*np.array( [1]*len( mesh.x ) )
u = spsolve( M, Source )

pp.tricontourf( mesh, u )
mesh.plot_boundaries( labels = False, color = 'black' )

##########################################################################################################

pp.xlim(x_lim); pp.ylim(y_lim)
pp.axis('off')
pp.xticks([])
pp.yticks([])

pp.savefig( '../figures/' + __file__.split('/')[-1].split('.')[0] + '_2.svg' , bbox_inches = 'tight' )

pp.show()

##########################################################################################################
