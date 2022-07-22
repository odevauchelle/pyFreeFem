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

for _ in range(4) :

    try :
        Th = pyff.adaptmesh( Th, v, iso = 1 )
    except :
        pass


    z = Th.x + 1j*Th.y
    v = real( sqrt( -1j*( z - 1j*y_wire[-1] ) ) ) + Th.x


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


script = pyff.InputScript( Th = Th )
script += pyff.edpScript('''
fespace Vh0( Th, P0 );
fespace Vh1( Th, P1 );

Vh0 area;
varf varea (unused, v) = int2d(Th)(v);
area[] = varea(0, Vh0);
''')

variational_forms = dict( dxv = 'int2d(Th)( u*dx(v) )', dyv = 'int2d(Th)( u*dy(v) )' )

script += pyff.VarfScript( **variational_forms, fespaces = ('Vh0', 'Vh1') )
script += pyff.OutputScript( area = 'vector' )

matrices = script.get_output()

dxv = v*matrices['dxv']/matrices['area']
dyv = v*matrices['dyv']/matrices['area']

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

name_to_index, index_to_name = Th.get_boundary_label_conversion()

from scipy.sparse import lil_matrix

flux = lil_matrix( shape( matrices['dxv'] )[::-1], dtype=np.cfloat )

node_indices = Th.get_boundaries()['wire'][0] # This is dangerous, as ordering can be messed up. We should keep track of the wire independently from Th.
node_indices = node_indices[::-1] # start from tip

for i in range( 1, len( node_indices ) ) :

    # identify the two triangles which share edge ( i-1, i )

    edge = [ node_indices[i-1], node_indices[i] ]
    dz = diff( Th.x[edge] + 1j*Th.y[edge] )[0]

    the_sign = -1

    for edge in ( edge, edge[::-1] ) : # two triangles per edge !
        triangle_index = pyff.find_triangle_index( Th.triangles, *edge )
        flux[ triangle_index, node_indices[i] ] = the_sign*dz/matrices['area'][triangle_index]
        the_sign *= -1

q_mat = matrices['dxv']*imag( flux ) - matrices['dyv']*real(flux)
# print(shape(flux))

q = v*q_mat
q = q[ node_indices ]
Q = cumsum(q)


##########################
#
# Plot discharge
#
##########################

figure()

ax_q = gca()

y = Th.y[ node_indices ]

y_th = linspace( min(y), max(y), 150 )
r = abs( y_th - y_th[-1] )

ax_q.plot( y, Q, '.')

ax_q.plot( y_th, 2*sqrt(r),'-')

show()
