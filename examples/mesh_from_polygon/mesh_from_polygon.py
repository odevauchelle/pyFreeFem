from pylab import *
import sys

sys.path.append('/home/olivier/git/pyFreeFem/')

import pyFreeFem as pyff

theta = linspace(0,2*pi,51)
r = 1 + .3*cos(5*theta)
x = r*cos(theta)
y = r*sin(theta)

plot(x,y,'o', color = 'grey')

# Th = pyff.TriMesh_from_polygon( array( [ x, y ] ).T )

labels = 'top', 'bottom'
i_cut = int( len(theta)/2 )
points = array([x,y]).T
boundaries = [ points[:i_cut], points[i_cut:]  ]

Th = pyff.TriMesh_from_boundaries( boundaries, labels )

Th.plot_triangles( color = 'lightgrey' )


axis('scaled')
axis('off')
xticks([])
yticks([])

Th.plot_boundaries()
legend()

# savefig( 'mesh_from_polygon_3.svg', bbox_inches = 'tight' )


show()
