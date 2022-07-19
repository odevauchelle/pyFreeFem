
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
# Th.add_boundary_edges( [ len( box_points ) ] + list( arange( len( box_points ) ) ) + [ len( box_points ) ] , 'box' )
Th.add_boundary_edges( range( len( box_points ) ) , 'box' )
Th.add_boundary_edges( range( len( box_points ), len( box_points ) + len( wire_points ) ) , 'wire' )


#########################
#
# Plot mesh
#
#########################

ax_mesh = gca()

Th.plot_triangles( color = 'grey', labels = 'index', ax = ax_mesh )
Th.plot_boundaries(ax = ax_mesh)
Th.plot_nodes( color = 'grey', labels = 'index' , ax = ax_mesh )

legend(loc = 'upper center')

ax_mesh.axis('scaled')
ax_mesh.axis('off')
ax_mesh.set_xticks([])
ax_mesh.set_yticks([])


# savefig( '../../figures/wire_mesh.svg', bbox_inches = 'tight')


#########################
#
# Matrices
#
#########################

script = pyff.InputScript( Th = Th )

script += pyff.edpScript('fespace Vh( Th, P1 );')

variational_forms = dict( stiffness = 'int2d(Th)( dx(u)*dx(v) +  dy(u)*dy(v) )' )

name_to_index, index_to_name = Th.get_boundary_label_conversion()

for name, index in name_to_index.items() :
    variational_forms.update( { 'BoundaryGramian_' + name : 'int1d(Th,' + str(index) + ')( u*v )' } )

print(variational_forms)

script += pyff.VarfScript( **variational_forms )

matrices = script.get_output()

#########################
#
# Look at matrices
#
#########################

fig, axs = subplots( ncols = len( matrices ), figsize = (8,3), sharey = True )

for i, ( name, mat ) in enumerate( matrices.items() ) :

    ax = axs[i]
    ax.imshow( mat.toarray() )
    ax.set_title(name)

# savefig( '../../figures/wire_matrices.svg', bbox_inches = 'tight')

#########################
#
# Absorbing boundary conditions
#
#########################

from scipy.sparse.linalg import spsolve

epsilon = 1e-6
ones_vector = Th.x*0 + 1.

v = spsolve(
    matrices['stiffness'] + 1/epsilon*( matrices['BoundaryGramian_box'] + matrices['BoundaryGramian_wire'] ),
    matrices['BoundaryGramian_box']*ones_vector
    )

figure()
ax_v = gca()

ax_v.tricontourf( Th, v )
ax_v.axis('scaled')
ax_v.axis('off')

Th.plot_boundaries(ax = ax_v, clip_on = False)
Th.plot_triangles(ax = ax_v, color = 'w', alpha = .3 )

ax_v.set_xticks([])
ax_v.set_yticks([])

# savefig( '../../figures/wire_field.svg', bbox_inches = 'tight')


show()
