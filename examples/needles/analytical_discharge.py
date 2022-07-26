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

y_wire = linspace(0,.5,3)

wire_points = list( array( [ y_wire*0, y_wire ] ).T )
# wire_points = wire_points[::-1][1:] + wire_points


Th = pyff.TriMesh( *array( box_points + wire_points ).T )
Th.add_boundary_edges( range( len( box_points ), len( box_points ) + len( wire_points ) ) , 'wire' )

fig_domain = figure()
ax_domain = gca()
Th.plot_triangles( ax = ax_domain, lw = .7, color = 'grey' )
Th.plot_boundaries( ax = ax_domain )

ax_domain.axis('scaled')
ax_domain.axis('off')
ax_domain.set_xticks([])
ax_domain.set_yticks([])
ax_domain.legend()

# savefig( '../../figures/' + __file__.split('/')[-1].split('.')[0] + '_mesh.svg' , bbox_inches = 'tight' )

#########################
#
# field
#
#########################

def field( Th, a = 0 ) :
    z = Th.x + 1j*Th.y
    return real( sqrt( -1j*( z - 1j*y_wire[-1] ) ) - a*z )


#########################
#
# refine the mesh
#
#########################

for _ in range(3) :

    try :
        Th = pyff.adaptmesh( Th, v, iso = 1 )
    except :
        pass

    v = field( Th )

#########################
#
# Plot field
#
#########################

a_list = [0, 1, 1j]

fig_v, axs = subplots( ncols = len( a_list ), figsize = array( [ len( a_list ), 1.1 ] )*3 )

for i, a in enumerate(a_list) :

    ax = axs[i]

    ax.tricontourf( Th, field( Th, a ) )

    Th.plot_triangles( ax = ax, lw = .7, color = 'w', alpha = .2 )
    Th.plot_boundaries( ax = ax )

    ax.axis('scaled')
    ax.axis('off')
    ax.set_xticks([])
    ax.set_yticks([])
    ax.set_title( r'$a=$' + str(a) )

legend(loc = 'upper center')

# savefig( '../../figures/' + __file__.split('/')[-1].split('.')[0] + '_field.svg' , bbox_inches = 'tight' )


########################
#
# Failure of int1d
#
########################

name_to_index, index_to_name = Th.get_boundary_label_conversion()

script = pyff.InputScript( Th = Th )
script += 'fespace Vh(Th,P1);'
script += pyff.InputScript( v = 'vector' )
script += 'real Q = int1d( Th,' + str( name_to_index['wire'] ) + ' )( N.x*dx(v) + N.y*dy(v) );'
script += pyff.OutputScript( Q = 'real' )

for i, a in enumerate(a_list) :
    print( 'a=', a, script.get_output( v = field( Th, a ) ) )


#########################
#
# Gradient
#
#########################

matrices = pyff.gradient_matrices( Th )

dxv = matrices['grad_x']*v
dyv = matrices['grad_y']*v

print( 'v', len(Th.x), len(v) )
print( 'dxv', len(Th.triangles), len(dxv) )

figure()
ax_grad = gca()

X = []

for i, triangle in enumerate( Th.triangles ) :

    X += [ [ mean( Th.x[triangle] ), mean( Th.y[triangle] ) ]  ]


ax = ax_grad
ax.tricontourf( Th, field( Th, a ) )
Th.plot_triangles( ax = ax, lw = .7, color = 'w', alpha = .2 )
Th.plot_boundaries( ax = ax )

ax.quiver( *array(X).T, -dxv, -dyv )

ax.axis('scaled')
ax.axis('off')
ax.set_xticks([])
ax.set_yticks([])
ax.set_title( r'$a=$' + str(a) )

# savefig( '../../figures/' + __file__.split('/')[-1].split('.')[0] + '_gradient.svg' , bbox_inches = 'tight' )

#######################
#
# Flux along wire
#
######################

boundary_nodes = Th.get_boundaries()['wire'][0] # This is dangerous, as ordering can be messed up. We should keep track of the wire independently from Th.
boundary_nodes = boundary_nodes[::-1] # start from tip

q_mat = pyff.flux_along_needle( Th, boundary_nodes )

print( shape(q_mat) )

figure()

ax_q = gca()

y = Th.y[boundary_nodes]

ms = 15

for a in a_list:
    qds = q_mat.dot( field( Th, a ) )
    q = -qds[boundary_nodes][1:]/diff(y)

    ax_q.plot( y[1:], q, '.', ms = ms, label = '$a=$' + str(a) )
    ms *= .7

y_th = linspace( 0 , .5, 100 )
ax_q.plot( y_th, 1/sqrt( .5 - y_th ), color = 'grey', zorder = -1, label = 'exact' )


ax_q.legend()
ax_q.set_xlabel('$y$')
ax_q.set_ylabel('Flux into wire $q$')

# savefig( '../../figures/' + __file__.split('/')[-1].split('.')[0] + '_q.svg' , bbox_inches = 'tight' )

####################################
#
# dicharge along wire
#
####################################

print('Q=', sum(qds) )

Q_mat = pyff.needle_discharge( Th, boundary_nodes )

##########################
#
# Plot discharge
#
##########################

figure()

ax_Q = gca()

ax_Q.plot( y, Q_mat.dot( field( Th, 0 ) )[boundary_nodes], '.')
ax_Q.plot( y_th, 2*sqrt(y_th[-1] - y_th), color = 'grey', label = 'exact' )

ax_Q.set_xlabel('$y$')
ax_Q.set_ylabel('Discharge $Q$')

ax_Q.legend()
# savefig( '../../figures/' + __file__.split('/')[-1].split('.')[0] + '_Q.svg' , bbox_inches = 'tight' )


show()
