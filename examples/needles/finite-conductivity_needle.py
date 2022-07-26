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
wire_points = [ [-.1,0],[ -.1,.6] ]

# plot( *array(box_points).T, 'o', color = 'grey')
# plot( *array(wire_points).T, 'o', color = 'grey')

Th = pyff.TriMesh( *array( box_points + wire_points ).T )
# Th.add_boundary_edges( [ len( box_points ) ] + list( arange( len( box_points ) ) ) + [ len( box_points ) ] , 'box' )
Th.add_boundary_edges( range( len( box_points ) )[1:-1] , 'box' )
Th.add_boundary_edges( range( len( box_points ), len( box_points ) + len( wire_points ) ) , 'wire' )
Th.add_boundary_edges( [3,4,0], 'bottom' )


#########################
#
# Matrices
#
#########################

script = pyff.InputScript( Th = 'mesh' )

script += pyff.edpScript('fespace Vh( Th, P1 );')

variational_forms = dict( stiffness = 'int2d(Th)( dx(u)*dx(v) +  dy(u)*dy(v) )' )

name_to_index, index_to_name = Th.get_boundary_label_conversion()

for name, index in name_to_index.items() :
    variational_forms.update( { 'BoundaryGramian_' + name : 'int1d(Th,' + str(index) + ')( u*v )' } )


variational_forms['wire_u_dv_ds'] = 'int1d(Th,' + str( name_to_index['wire'] ) + ')( v*N.x*dy(u) - v*N.y*dx(u) )'
# variational_forms['wire_u_dv_ds'] = 'int1d(Th,' + str( name_to_index['wire'] ) + ')( u*N.x*dy(v) - u*N.y*dx(v) )'

script += pyff.VarfScript( **variational_forms )

#########################
#
# Absorbing boundary conditions
#
#########################

from scipy.sparse.linalg import spsolve

epsilon = 1e-6
kappa = 2

for _ in range(4) :

    try :
        Th = pyff.adaptmesh( Th, v, err = 3e-3 )
    except :
        pass

    matrices = script.get_output( Th = Th )


    ones_vector = Th.x*0 + 1.


    v = spsolve(
        matrices['stiffness'] + 1/epsilon*( matrices['BoundaryGramian_box'] + matrices['BoundaryGramian_wire'] + matrices['BoundaryGramian_bottom']),
        1/epsilon*matrices['BoundaryGramian_box']*ones_vector
        )

boundary_nodes = Th.get_boundaries()['wire'][0] # This is dangerous, as ordering can be messed up. We should keep track of the wire independently from Th.
boundary_nodes = array(boundary_nodes)[ argsort( Th.y[ boundary_nodes ] ) ][::-1]

Q_mat = pyff.needle_discharge( Th, boundary_nodes )


v = spsolve(
    matrices['stiffness'] + 1/epsilon*( matrices['BoundaryGramian_box'] + matrices['BoundaryGramian_bottom'] + matrices['BoundaryGramian_wire'].dot(Q_mat) - kappa*matrices['wire_u_dv_ds'] ),
    1/epsilon*matrices['BoundaryGramian_box']*ones_vector
    )

figure()
ax_v = gca()

v_contours = ax_v.tricontourf( Th, v )
# colorbar(v_contours)

ax_v.axis('scaled')
ax_v.axis('off')

Th.plot_boundaries(ax = ax_v, clip_on = False)
# Th.plot_triangles(ax = ax_v, color = 'w', alpha = .3, lw = 1 )

ax_v.set_xticks([])
ax_v.set_yticks([])

show()
