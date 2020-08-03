from pylab import *

import sys
sys.path.append('../../')
import pyFreeFem as pyff

#############################
#
# GENERIC SCRIPT
#
#############################

script = pyff.InputScript( Th = 'mesh' )

script += '''
// Fespace

fespace Vh(Th, P1);
Vh phi, w;

// Solve
problem Poisson(phi, w)
    = int2d(Th)(
          dx(phi)*dx(w)
        + dy(phi)*dy(w)
    )
    - int2d(Th)(
          w
    )
    + on(1, phi=0)
    ;

for( int i=1; i<4; i++ )
    {
    Poisson;
    Th = adaptmesh( Th, phi, iso = 1 ) ;
    }

Poisson;

real discharge;
discharge = int2d(Th)(phi);
'''

script += pyff.OutputScript( phi = 'vector', Th = 'mesh', discharge = 'real' )


##############################
#
# Function
#
##############################

def get_discharge( W, D = 1 ) :

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

    FF_out = script.get_output( Th = Th )

    return FF_out['Th'], FF_out['phi'], FF_out['discharge']

Th, u, discharge = get_discharge( W = .2 )

tricontourf(Th,u)
Th.plot_triangles(lw = .5, color = 'w')
Th.plot_boundaries( zorder = 10, clip_on = False)

ax = gca()
ax.axis('equal'); ax.axis('off')
ax.set_xticks([]); ax.set_yticks([])

savefig( '../../figures/' + __file__.split('/')[-1].split('.')[0] + '.svg' , bbox_inches = 'tight' )


show()
