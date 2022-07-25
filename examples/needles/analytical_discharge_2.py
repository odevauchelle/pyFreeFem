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

y_wire = linspace(0,.5,30)

wire_points = list( array( [ y_wire*0, y_wire ] ).T )
# wire_points = wire_points[::-1][1:] + wire_points


Th = pyff.TriMesh( *array( box_points + wire_points ).T )
Th.add_boundary_edges( range( len( box_points ), len( box_points ) + len( wire_points ) ) , 'wire' )

for _ in range(10) :

    try :
        Th = pyff.adaptmesh( Th, v, iso = 1 )
    except :
        pass


    z = Th.x + 1j*Th.y
    v = real( sqrt( -1j*( z - 1j*y_wire[-1] ) ) ) - Th.x


#########################
#
# Plot field
#
#########################

ax_mesh = gca()

tricontourf( Th, v )

Th.plot_triangles( ax = ax_mesh, lw = 1, color = 'w', alpha = .2 )
Th.plot_boundaries( ax = ax_mesh )

legend(loc = 'upper center')

ax_mesh.axis('scaled')
ax_mesh.axis('off')
ax_mesh.set_xticks([])
ax_mesh.set_yticks([])

#########################
#
# Gradient
#
#########################

from scipy.sparse.linalg import spsolve

matrices = pyff.gradient_matrices( Th )

dxv = matrices['grad_x']*v
dyv = matrices['grad_y']*v

x = []
y = []

for i, triangle in enumerate( Th.triangles ) :

    x += [ mean( Th.x[triangle] )]
    y += [ mean( Th.y[triangle] )]

ax_mesh.quiver( x, y, -dxv, -dyv )

####################################
#
# flux along wire
#
####################################

boundary_nodes = Th.get_boundaries()['wire'][0] # This is dangerous, as ordering can be messed up. We should keep track of the wire independently from Th.
boundary_nodes = boundary_nodes[::-1] # start from tip

Q_mat = pyff.needle_discharge( Th, boundary_nodes )
Q = Q_mat*v

##########################
#
# Plot discharge
#
##########################

figure()

ax_q = gca()

y = Th.y[ boundary_nodes ]

y_th = linspace( min(y), max(y), 150 )
r = abs( y_th - y_th[-1] )

ax_q.plot( y, Q[boundary_nodes], '.')

ax_q.plot( y_th, 2*sqrt(r),'-')

show()
