from pylab import *

import sys
sys.path.append('./../../')

import pyFreeFem as pyff


theta =  linspace( 0, pi, 25 )
top = cos(theta), 0.5*sin(theta)

theta =  linspace( -pi, 0, len(theta) )[1:-1]
bottom = cos(theta), sin(theta)

x = list( top[0] ) + list( bottom[0] ) + [0]
y = list( top[1] ) + list( bottom[1] ) + [0]

triangles = [ [ i, i + 1, len(x) - 1 ] for i in range( len(x) - 2 ) ]
triangles += [ [ len(x) - 2, 0, len(x) - 1 ] ]

Th = pyff.TriMesh( x, y, triangles )
Th.add_boundary_edges( range( len( top[0] ) ), 'top' )
Th.add_boundary_edges( list( range( len( top[0] ) - 1, len( x ) - 1 ) ) + [0], 'bottom' )

Th = pyff.adaptmesh( Th )

Th.plot_triangles( lw = 1, alpha = .5 )
Th.plot_boundaries()

legend()
axis('equal')
axis('off')
xticks([]), yticks([])
show()
