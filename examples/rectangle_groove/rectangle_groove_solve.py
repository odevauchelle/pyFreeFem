from pylab import *

import sys
sys.path.append('../../')
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
# add boundaries
#
##############################

Th.add_boundary_edges( [1,2,3,0] )
Th.add_boundary_edges( [0,1] )

##############################
#
# Refine mesh
#
##############################

Th = pyff.adaptmesh(Th, hmax = .03, iso = 1)

##############################
#
# Solve
#
##############################

script = pyff.InputScript( Th = Th )

print(script.get_edp())

script += '''
// Fespace
fespace Vh(Th, P1);
Vh phi, w;

// Solve
solve Poisson(phi, w)
    = int2d(Th)(
          dx(phi)*dx(w)
        + dy(phi)*dy(w)
    )
    - int2d(Th)(
          w
    )
    + on(1, phi=0)
    ;
'''

script += pyff.OutputScript( phi = 'vector' )

u = script.get_output()['phi']

tricontourf(Th,u)
Th.plot_boundaries( zorder = 10, clip_on = False)

ax = gca()
ax.axis('equal'); ax.axis('off')
ax.set_xticks([]); ax.set_yticks([])

savefig( '../../figures/' + __file__.split('/')[-1].split('.')[0] + '.svg' , bbox_inches = 'tight' )


show()
