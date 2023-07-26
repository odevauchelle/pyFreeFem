
from pylab import *
from scipy.sparse.linalg import svds

import sys
sys.path.append('./../../')

import pyFreeFem as pyff

#################################
#
# build mesh
#
################################

h = .05

theta = linspace( 0, 2*pi, int(2*pi/h))[:-1]

external_boundary = array([ cos(theta), sin(theta) ]).T

r = 0.3
X_b = array( [ -0.5, 0 ] )

theta = linspace( 0, 2*pi, int( r*2*pi/h ) )[:-1]
internal_boundary =  array([ cos(theta), sin(theta) ])[::-1].T*r + X_b

plot( *external_boundary.T, 'o' )
plot( *internal_boundary.T, 'o' )

xy = list( external_boundary )
external_edges = [ [ i, i + 1 ] for i in range( len( external_boundary ) - 1 ) ]
external_edges += [ [ len(external_boundary) - 1 , 0 ] ]

xy += list( internal_boundary )
internal_edges = [ [ i, i + 1 ] for i in range( len( internal_boundary ) - 1 ) ]
internal_edges += [ [ len(internal_boundary) - 1 , 0 ] ]
internal_edges = list( array( internal_edges ) + len(external_boundary) )

Th = pyff.TriMesh( *array(xy).T )
Th.add_boundary_edges( external_edges, 'external' )
Th.add_boundary_edges( internal_edges, 'internal' )

Th = pyff.adaptmesh( Th, hmax = h )

Th.plot_triangles( color = 'k', alpha = .3, lw = 1 )
Th.plot_boundaries()
axis('equal'); axis('off')

#################################
#
# temperature field & matrices
#
################################

region = array( [0]*len( Th.triangles ) )
triangle_centers = []


for i in range( len( Th.triangles ) ) :

    xy_center = mean( Th.x[ Th.triangles[i] ] ), mean( Th.y[ Th.triangles[i] ] )

    triangle_centers += [ xy_center ]

    if norm( xy_center - X_b ) < r :
        region[i] = 1
        plot( *xy_center, '.', color = 'grey' )

script = pyff.InputScript( Th = Th )

script +='''
fespace Vh0( Th, P0 );
fespace Vh( Th, P1 );
Vh0 region;

Vh u, v;

'''

script += pyff.InputScript( region = region, declare = False )

script += pyff.VarfScript(
    stiffness = 'int2d(Th)( dx(u)*dx(v) +  dy(u)*dy(v) )',
    gravity_in = 'int2d(Th)( dy(u)*v*region )',
    gravity_out = 'int2d(Th)( dy(u)*v*( 1 - region ) )',
    )

script += pyff.VarfScript( fespaces = ('Vh','Vh0'),
    gramianP0 = 'int2d(Th)( u*v )',
    grad_x = 'int2d(Th)( dx(u)*v )',
    grad_y = 'int2d(Th)( dy(u)*v )',
    )

# script += '''
# plot(region,fill=1,ps="region.eps");
# '''

ff_out = script.get_output()

#################################
#
# solve
#
################################

beta_out = 1
beta_in = beta_out*1.2
beta = beta_in*region + beta_out*( 1 - region )
force = 1

lsv, v, rsv = svds( ff_out['stiffness'] + force*( ff_out['gravity_in']*beta_in + ff_out['gravity_out']*beta_out ) , k = 1, which = 'SM')
print(v)
rho_D = lsv[:,0]

################################
#
# normalization
#
################################

D_out = 1
D_in = D_out*( beta_in/beta_out )**(-3/2)
D = D_in*region + D_out*( 1 - region )

Z = 1/D.T@( ff_out['gramianP0']@rho_D )

rho_D /= Z

proj01 = pyff.get_projector( Th, 'P0', 'P1' )

rho = rho_D/( proj01@D )

################################
#
# current
#
################################
nabla = pyff.gradient_matrices( Th, 'P0', 'P1' )
proj10 = pyff.get_projector( Th, 'P1', 'P0' )

qx = -nabla['grad_x']@rho_D
qy = -nabla['grad_y']@rho_D - force*beta*( proj10@rho_D )


# figure()
#
# plot( array( triangle_centers )[:,1], force*beta*( proj10@rho_D ), '.' )
# plot( array( triangle_centers )[:,1], -ff_out['grad_y']@rho_D, '.' )

figure()

tricontourf( Th, rho_D )
Th.plot_boundaries()

quiver( *array( triangle_centers ).T, qx, qy  )

axis('equal'); axis('off')

show()
