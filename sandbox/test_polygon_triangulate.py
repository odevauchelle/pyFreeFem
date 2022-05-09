from pylab import *
from matplotlib.tri import Triangulation

from polygon_triangulate import polygon_triangulate

theta = linspace(0,2*pi,51)[:-1]
r = 1 + .3*cos(5*theta)
x = r*cos(theta)
y = r*sin(theta)

triangles = polygon_triangulate( len(x), x, y )
# triangles[triangles==-1] = 0
# print(triangles)
tri = Triangulation( x, y, triangles = triangles )
#
triplot(tri, color = 'lightgrey')
plot( x, y, '--o' )
axis('scaled')
axis('off')
show()
