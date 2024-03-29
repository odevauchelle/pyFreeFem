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
wire_points = [ [-.1,0], [-.1,.6] ]

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

script += pyff.VarfScript( **variational_forms )

matrices = script.get_output( Th = Th )

#########################
#
# Look at matrices
#
#########################

fig, axs = subplots( ncols = len( matrices ), figsize = array([len(matrices),1])*3, sharey = True )

for i, ( name, mat ) in enumerate( matrices.items() ) :

    ax = axs[i]
    ax.imshow( mat.toarray() )
    ax.set_title(name)

savefig( '../../figures/wire_matrices.svg', bbox_inches = 'tight')

#########################
#
# Absorbing boundary conditions
#
#########################

from scipy.sparse.linalg import spsolve

epsilon = 1e-6

for _ in range(4) :

    try :
        Th = pyff.adaptmesh( Th, v, err = 3e-3 )
        matrices = script.get_output( Th = Th )
    except :
        pass

    ones_vector = Th.x*0 + 1.

    v = spsolve(
        matrices['stiffness'] + 1/epsilon*( matrices['BoundaryGramian_box'] + matrices['BoundaryGramian_wire'] + matrices['BoundaryGramian_bottom']),
        1/epsilon*matrices['BoundaryGramian_box']*ones_vector
        )

#########################
#
# Plot absorbing field
#
#########################

figure()
ax_v = gca()

v_contours = ax_v.tricontourf( Th, v )
# colorbar(v_contours)

ax_v.axis('scaled')
ax_v.axis('off')

Th.plot_boundaries(ax = ax_v, clip_on = False)
Th.plot_triangles(ax = ax_v, color = 'w', alpha = .2, lw = 0.7 )

ax_v.set_xticks([])
ax_v.set_yticks([])

# savefig( '../../figures/wire_field.svg', bbox_inches = 'tight')


###########################
#
# finite resisitivity
#
###########################

boundary_nodes = Th.get_boundaries()['wire'][0] # This is dangerous, as ordering can be messed up. We should keep track of the wire independently from Th.
boundary_nodes = array(boundary_nodes)[ argsort( Th.y[ boundary_nodes ] ) ][::-1]

Q_mat = pyff.needle_discharge( Th, boundary_nodes )


script += pyff.VarfScript( wire_u_dv_ds = 'int1d(Th,' + str( name_to_index['wire'] ) + ')( v*N.x*dy(u) - v*N.y*dx(u) )' )
matrices = script.get_output( Th = Th )


kappa_list = logspace( -1, 1, 3  )

fig_fc, axs = subplots( ncols = len(kappa_list), figsize = array([len(kappa_list),1])*3 )

fig_Q = figure()
ax_Q = gca()

for i, kappa in enumerate( kappa_list ) :

    ax = axs[i]

    v = spsolve(
        matrices['stiffness'] + 1/epsilon*( matrices['BoundaryGramian_box'] + matrices['BoundaryGramian_bottom'] + matrices['BoundaryGramian_wire'].dot(Q_mat) - kappa*matrices['wire_u_dv_ds'] ),
        1/epsilon*matrices['BoundaryGramian_box']*ones_vector
        )

    ax.tricontourf( Th, v )

    ax.axis('scaled')
    ax.axis('off')

    label = r'$\kappa=$' + str(kappa)

    ax.set_title( label )
    Th.plot_boundaries(ax = ax, clip_on = False)

    ax.set_xticks([])
    ax.set_yticks([])

    y = Th.y[boundary_nodes]
    Q = Q_mat.dot(v)[boundary_nodes]

    ax_Q.plot( y, Q, label = label )

ax_Q.legend()
ax_Q.set_xlabel('$y$')
ax_Q.set_ylabel('Discharge $Q$')

# fig_fc.savefig( '../../figures/wire_kappa.svg', bbox_inches = 'tight')
# fig_Q.savefig( '../../figures/wire_Q.svg', bbox_inches = 'tight')

show()
