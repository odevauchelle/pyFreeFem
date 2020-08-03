from pylab import *

import sys
sys.path.append('../')
import pyFreeFem as pyff

##############################
#
# PARAMETERS
#
##############################

D = 0.3
W = 1

##############################
#
# MESH
#
##############################

x = array([ .5,-.5,-.5,.5 ])*W
y = array([0,0,-1,-1])*D

Th = pyff.TriMesh( x, y )

##############################
#
# initial plot
#
##############################

figure()
ax = gca()

Th.plot_triangles(ax = ax)
Th.plot_nodes( labels = 'index', ax = ax, color = 'tab:blue' )

ax.axis('equal'); ax.axis('off')
ax.set_xticks([]); ax.set_yticks([])

savefig( '../figures/' + __file__.split('/')[-1].split('.')[0] + '_initial_mesh.svg' , bbox_inches = 'tight' )

##############################
#
# add boundaries
#
##############################
Th.add_boundary_edges( [1,2,3,0] )
Th.add_boundary_edges( [0,1] )

figure()
ax = gca()

Th.plot_triangles(ax = ax)
Th.plot_boundaries(ax = ax)
ax.legend(title = 'Boundaries')

Th.plot_nodes( labels = 'index', ax = ax, color = 'tab:blue' )

ax.axis('equal'); ax.axis('off')
ax.set_xticks([]); ax.set_yticks([])

##############################
#
# Refine mesh
#
##############################

Th = pyff.adaptmesh(Th, hmax = .03, iso = 1)

figure()
ax = gca()

Th.plot_triangles(ax = ax)
Th.plot_boundaries(ax = ax)
ax.legend(title = 'Boundaries')

ax.axis('equal'); ax.axis('off')
ax.set_xticks([]); ax.set_yticks([])



show()

# script = pyff.InputScript( D = D, W = W )
