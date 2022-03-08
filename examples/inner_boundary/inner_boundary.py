# coding:utf8
from pylab import *
from matplotlib.tri import Triangulation

import sys
sys.path.append('../../')
import pyFreeFem as pyff

# geometries

# box
Xb = array([[0.5,0,0,1,1], [0,0,1,1,0] ]).T

# river
yr = linspace(0,0.7,30) - .0
xr = 0.5 + 0.1*cos(10*yr)
Xr = array([xr,yr]).T

figure()
plot( *Xb.T,'o', color='grey')
plot( *Xr.T,'.', color='tab:blue')

axis('scaled')
axis('off')
xticks([]); yticks([])

# show()

#topologies

divide = [ [i,i+1] for i in range(len(Xb)-1) ]
divide += [ [ len(Xb) - 1, 0  ] ]

river = [ [i,i+1] for i in range(len(Xr)-1) ]

# merge river and box

# 1 geometric merging

river_mouth_index = 0 # node
divide_edge_index = len(divide) - 1 # edge

X = array( list(Xr) + list(Xb) )
divide = ( array(divide) + len(Xr) ).tolist()

# topologic merging

removed_edge = divide.pop( divide_edge_index )
new_edges = [ [ removed_edge[0], river_mouth_index ], [river_mouth_index, removed_edge[1]] ]

for new_edge in new_edges[::-1] :
    divide.insert( divide_edge_index, new_edge )

# create initial triangulation

tri = Triangulation(*X.T)
triangles = tri.triangles

# print(triangles,len(triangles),len(x))
Th = pyff.TriMesh( *X.T, triangles)
# Th.add_boundary_edges(  )

Th.add_boundary_edges( river, label = 'river' )
Th.add_boundary_edges( divide, label = 'divide' )

Th = pyff.adaptmesh(Th, hmax = .05, iso = 1)
#

Th.plot_triangles(color = 'k', alpha = .3)

Th.plot_boundaries()
legend()

show()
