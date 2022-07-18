
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

#########################
#
# Plot mesh
#
#########################


Th.plot_triangles( color = 'grey', labels = 'index' )
Th.plot_boundaries()
Th.plot_nodes( color = 'grey', labels = 'index' )

legend(loc = 'upper center')

axis('scaled')
axis('off')


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

show()
