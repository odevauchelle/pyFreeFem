import matplotlib.pyplot as pp
from tempfile import NamedTemporaryFile
from random import sample

import sys
sys.path.append('./../')

import pyFreeFem as pyff

# Create mesh with FreeFem++

script = pyff.edpScript( '''
border Circle( t = 0, 2*pi ){ x = cos(t); y = sin(t); }
mesh Th = buildmesh( Circle(150) );
''' )


script += pyff.edpOutput( type = 'mesh', name = 'Th' )

Th = script.get_output()['Th']


# Change mesh
for triangle_index in sample( range( len( Th.triangles ) ), 20 ) :
    Th.boundary_edges.update( { ( triangle_index, 0 ) : 2 } )

Th.rename_boundary( {1:'initial', 2:'new'} )

# Export mesh back to FreeFem

script = pyff.edpScript( pyff.edpInput( name = 'Th', source = Th ) )

# calculate FEM matrices

script +='''
fespace Vh( Th, P1 ) ;
Vh u, v ;
'''

matrices = {
    'stiffness' : 'int2d(Th)( dx(u)*dx(v) +  dy(u)*dy(v) )',
    'Grammian' : 'int2d(Th)( u*v )',
    'boundary_Grammian' : 'int1d(Th, 1, 2)( u*v )'
}

for matrix_name in matrices.keys() :
    script += pyff.VarfBlock( name = matrix_name, varf = matrices[matrix_name] )

matrices = script.get_output()

##########################################################################################################

from scipy.sparse.linalg import spsolve
import numpy as np

epsilon = 1e-4
M = - matrices[ 'stiffness' ] + 1./epsilon*matrices[ 'boundary_Grammian' ]
Source = matrices[ 'Grammian' ]*np.array( [1]*len( Th.x ) )
u = spsolve( M, Source )


##################
#
# First FIGURE
#
#################

figs = {}

figs['mesh'] = pp.figure()
Th.plot_triangles( color = 'k', lw = .5, alpha = .2  )
Th.plot_boundaries()
pp.legend(title = 'Boundary')
x_lim, y_lim = pp.gca().get_xlim(), pp.gca().get_ylim()

figs['field'] = pp.figure()
pp.tricontourf( Th, u )
Th.plot_boundaries( color = 'black' )


for fig_key in figs.keys() :

    pp.figure(figs[fig_key].number)

    pp.axis('equal')
    pp.axis('off')
    pp.xticks([])
    pp.yticks([])

    pp.xlim(x_lim); pp.ylim(y_lim)

    pp.savefig( '../figures/' + __file__.split('/')[-1].split('.')[0] + '_' + fig_key + '.svg' , bbox_inches = 'tight' )

pp.show()
