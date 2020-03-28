from pylab import *
from scipy.sparse.linalg import spsolve

import sys
sys.path.append('./../../')

import pyFreeFem as pyff

###########################
#
#  PARAMETERS
#
############################

width = 5.
height = 1.
amplitude = .3
k = 3*2*pi/width
npts = 10

###########################
#
#  MESH
#
############################

script = pyff.InputScript( width = width, height = height, amplitude = amplitude, k = k, npts = npts )

script += pyff.edpScript('''
border top( t = width/2, -width/2 ){ x = t; y = height; }
border bottom( t = -width/2, width/2 ){ x = t; y = amplitude*cos( k*t ); }
border left( t = height,  amplitude*cos( -k*width/2 ) ){ x = -width/2; y = t; }
border right( t = amplitude*cos( k*width/2 ), height ){ x = width/2; y = t; }
mesh Th = buildmesh( top(npts) + bottom(10*npts) + left(npts) + right(npts) );
''')

script += pyff.OutputScript( Th = 'mesh' )

Th = script.get_output()['Th']

Th = pyff.adaptmesh( Th, hmax = height/15 )

###########################
#
#  FE MATRICES
#
############################

script = pyff.InputScript( Th = Th )

script += '''
fespace Vh1( Th, P1 );
fespace Vh2( Th, P2 );
'''

script += pyff.VarfScript( fespaces = ( 'Vh2', 'Vh2' ),
    stiffness = 'int2d(Th)( dx(u)*dx(v) +  dy(u)*dy(v) )',
    Gramian = 'int2d(Th)( u*v )',
    bottom_Gramian = 'int1d( Th, 2 )( u*v )'
    )

script += pyff.VarfScript( fespaces = ( 'Vh1', 'Vh1' ),
        Gramian_P1 = 'int2d(Th)( u*v )'
        )

script += pyff.VarfScript( fespaces = ( 'Vh2', 'Vh1' ),
        Dx = 'int2d(Th)( v*dx(u) )',
        Dy = 'int2d(Th)( v*dy(u) )',
        )

M = script.get_output()

###########################
#
#  SOLVE
#
############################

epsilon = 1e-6
Poisson = - M['stiffness'] + 1./epsilon*M['bottom_Gramian']
Poisson_source = M['Gramian']*np.array( [-1]*shape( M['Gramian'] )[0] )

u = spsolve( Poisson, Poisson_source )
dx_u = spsolve( M['Gramian_P1'], M['Dx']*u )
dy_u = spsolve( M['Gramian_P1'], M['Dy']*u )
tau = sqrt( dx_u**2 + dy_u**2 )

###########################
#
#  EXTRACT BOTTOM SHEAR STRESS
#
############################

bottom_label = 2

boundary_nodes = Th.get_boundaries()

###########################
#
#  PLOT
#
############################

tau_color = 'tab:red'
boundary_color = 'black'

fig, (ax_flow, ax_tau) = subplots( nrows = 2, sharex = True, figsize = 6*array([1,.6]) )

ax_flow.tricontourf( Th, pyff.interpolate( Th, u, 'P2', 'P1' ) )
Th.plot_triangles( ax = ax_flow, color = 'k', lw = .5, alpha = .2 )

for label in 1,3,4 :
    for nodes in  boundary_nodes[label] :
        ax_flow.plot( Th.x[ nodes ], Th.y[nodes], color = boundary_color, linestyle = '--', clip_on = False )

for nodes in  boundary_nodes[bottom_label] :
    ax_flow.plot( Th.x[ nodes ], Th.y[ nodes ], color = boundary_color, linestyle = '-', clip_on = False )
    ax_tau.plot( Th.x[nodes ], tau[nodes], color = tau_color )

ax_flow.axis('equal'); ax_flow.axis('off')
ax_flow.set_yticks([])

ax_tau.set_xlabel('Cross-stream coordinate $x$')
ax_tau.set_ylabel(r'Shear stress $\tau$')

# fig_path_and_name = './../../figures/' + __file__.split('/')[-1].split('.')[0] + '.svg'
# savefig( fig_path_and_name , bbox_inches = 'tight' )
# print(fig_path_and_name)

show()
