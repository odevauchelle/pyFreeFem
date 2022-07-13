
import sys
sys.path.append('./../../')

from pylab import *
import pyFreeFem as pyff

#########################
#
# Create mesh points
#
#########################

box_points = [ [ .5, 0 ], [.5,1], [-.5,1], [-.5,0] ]
wire_points = [ [0,0],[0,.2],[.1,.4] ]

# plot( *array(box_points).T, 'o', color = 'grey')
# plot( *array(wire_points).T, 'o', color = 'grey')

Th = pyff.TriMesh( *array( box_points + wire_points ).T )
Th.add_boundary_edges( [ len( box_points ) ] + list( arange( len( box_points ) ) ) + [ len( box_points ) ] , 'box' )
Th.add_boundary_edges( range( len( box_points ), len( box_points ) + len( wire_points ) ) , 'wire' )

Th.plot_triangles( color = 'grey', labels = 'index' )
Th.plot_boundaries()
Th.plot_nodes( color = 'grey', labels = 'index' )

legend(loc = 'upper center')

axis('scaled')
axis('off')

# savefig( '../../figures/wire_mesh.svg', bbox_inches = 'tight')

show()
