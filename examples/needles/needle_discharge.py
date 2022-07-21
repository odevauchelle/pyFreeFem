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
wire_points = [ [0,0],[0,.3],[.1,.6],[.15,.7] ]

# plot( *array(box_points).T, 'o', color = 'grey')
# plot( *array(wire_points).T, 'o', color = 'grey')

Th = pyff.TriMesh( *array( box_points + wire_points ).T )
# Th.add_boundary_edges( [ len( box_points ) ] + list( arange( len( box_points ) ) ) + [ len( box_points ) ] , 'box' )
Th.add_boundary_edges( range( len( box_points ) )[1:-1] , 'box' )
Th.add_boundary_edges( range( len( box_points ), len( box_points ) + len( wire_points ) ) , 'wire' )
Th.add_boundary_edges( [3,4,0], 'bottom' )


#########################
#
# Plot mesh
#
#########################

ax_mesh = gca()

Th.plot_triangles( color = 'grey', labels = 'index', ax = ax_mesh )
Th.plot_boundaries( ax = ax_mesh )
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

script = pyff.InputScript( Th = 'mesh' )

script += pyff.edpScript('fespace Vh( Th, P1 );')

variational_forms = dict( stiffness = 'int2d(Th)( dx(u)*dx(v) +  dy(u)*dy(v) )' )

name_to_index, index_to_name = Th.get_boundary_label_conversion()

for name, index in name_to_index.items() :
    variational_forms.update( { 'BoundaryGramian_' + name : 'int1d(Th,' + str(index) + ')( u*v )' } )

# print(variational_forms)
# variational_forms['wire_u_dv_ds'] = 'int1d(Th,' + str( name_to_index['wire'] ) + ')( - v*N.x*dy(u) + v*N.y*dx(u) )'
variational_forms['wire_u_dv_ds'] = 'int1d(Th,' + str( name_to_index['wire'] ) + ')( u*N.x*dy(v) - u*N.y*dx(v) )'

script += pyff.VarfScript( **variational_forms )

matrices = script.get_output( Th = Th )

#########################
#
# Look at matrices
#
#########################

fig, axs = subplots( ncols = len( matrices ), figsize = (10,3), sharey = True )

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

for _ in range(1) :

    try :
        Th = pyff.adaptmesh( Th, v, iso = 1 )
        matrices = script.get_output( Th = Th )
    except :
        pass

    ones_vector = Th.x*0 + 1.

    v = spsolve(
        matrices['stiffness'] + 1/epsilon*( matrices['BoundaryGramian_box'] + matrices['BoundaryGramian_wire'] + matrices['BoundaryGramian_bottom']),
        1/epsilon*matrices['BoundaryGramian_box']*ones_vector
        )



figure()
ax_v = gca()

v_contours = ax_v.tricontourf( Th, v )
# colorbar(v_contours)

ax_v.axis('scaled')
ax_v.axis('off')

Th.plot_boundaries(ax = ax_v, clip_on = False)
Th.plot_triangles(ax = ax_v, color = 'w', alpha = .3, lw = 1 )

ax_v.set_xticks([])
ax_v.set_yticks([])

# savefig( '../../figures/wire_field.svg', bbox_inches = 'tight')

#########################
#
# Current in wire
#
#########################

# we first need to build a P0 x P1 stiffness matrix

script = pyff.InputScript( Th = Th )
script += pyff.edpScript('''
fespace Vh0( Th, P0 );
fespace Vh1( Th, P1 );
''')

variational_forms = dict( h_dn_v = 'int1d(Th,' + str( name_to_index['wire'] ) + ')( u*N.x*dx(v) + u*N.y*dy(v) )' )

script += pyff.VarfScript( **variational_forms, fespaces = ('Vh0', 'Vh1') ) # u needs to be Vh0, and v to Vh1

matrices.update( script.get_output() )

# we then build a P1 x P0 Heaviside matrix

from scipy.sparse import lil_matrix

wire_Heaviside = lil_matrix( shape( matrices['h_dn_v'] )[::-1] )


print(shape(matrices['h_dn_v']))
nodes_indices, tri_indices = matrices['h_dn_v'].nonzero()

print('active triangles', set(tri_indices))
print('active nodes', set(nodes_indices))


node_indices = Th.get_boundaries()['wire'][0] # This is dangerous, as ordering can be messed up. We should keep track of the wire independently from Th.
node_indices = node_indices[::-1] # start from tip

triangle_indices = []

for i in range( 1, len( node_indices ) ) :

    # identify the two triangles which share edge ( i-1, i )

    edge = node_indices[i-1], node_indices[i]

    for edge in ( edge, edge[::-1] ) :
        triangle_indices += [ pyff.find_triangle_index( Th.triangles, *edge ) ]

    wire_Heaviside[ triangle_indices, node_indices[i] ] = 1

# print(wire_Heaviside.toarray())

Q_mat = - matrices['h_dn_v']*wire_Heaviside
Q_mat = Q_mat.transpose()

Q = Q_mat*v
# ax_v.tricontourf( Th, Q )
# print( Q[node_indices] )

#########################
#
# Finite-conductivity wire
#
#########################

kappa = 100

v = spsolve(
    matrices['stiffness'] + 1/epsilon*( matrices['BoundaryGramian_box'] + matrices['BoundaryGramian_bottom'] + Q_mat - kappa*matrices['wire_u_dv_ds'] ),
    1/epsilon*matrices['BoundaryGramian_box']*ones_vector
    )

ax_v.tricontourf( Th, v )


show()
