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

Th = pyff.TriMesh( *array( box_points + wire_points ).T )
Th.add_boundary_edges( range( len( box_points ), len( box_points ) + len( wire_points ) ) , 'wire' )

for _ in range(5) :

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
# Current in wire
#
#########################

# we first need to build a P0 x P1 stiffness matrix

name_to_index, index_to_name = Th.get_boundary_label_conversion()


script = pyff.InputScript( Th = Th )
script += pyff.edpScript('''
fespace Vh0( Th, P0 );
fespace Vh1( Th, P1 );
''')

variational_forms = dict( h_dn_v = 'int1d(Th,' + str( name_to_index['wire'] ) + ')( u*N.x*dx(v) + u*N.y*dy(v) )' )

script += pyff.VarfScript( **variational_forms, fespaces = ('Vh0', 'Vh1') ) # u needs to be Vh0, and v to Vh1

matrices = script.get_output()

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

ax_q.plot( y, Q[ node_indices ], '.')

ax_q.plot( y_th, sqrt(r),'-')


show()
