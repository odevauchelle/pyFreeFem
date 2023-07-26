
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

h = .15

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

for i in range( len( Th.triangles ) ) :

    xy_center = mean( Th.x[ Th.triangles[i] ] ), mean( Th.y[ Th.triangles[i] ] )

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

# script += '''
# plot(region,fill=1);
# '''

ff_out = script.get_output()

#################################
#
# solve (Boltzmann equilibrium)
#
################################

beta_out = 1
beta_in = beta_out*2

lsv, v, rsv = svds( ff_out['stiffness'] + ff_out['gravity_in']*beta_in + ff_out['gravity_out']*beta_out , k = 3, which = 'SM')
print(v)
rho_D = lsv[:,0]

figure()

tricontourf( Th, rho_D )
Th.plot_boundaries()

axis('equal'); axis('off')

show()
