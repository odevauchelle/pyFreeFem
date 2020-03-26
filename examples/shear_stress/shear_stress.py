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
amplitude = .2
k = 3*2*pi/width
npts = 3

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

Th = pyff.adaptmesh( Th, hmax = height/5 )

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

print( Th.boundary_edges )


###########################
#
#  PLOT
#
############################

figure()
ax = gca()

# u_1 = pyff.interpolate( Th, u, 'P2', 'P1' )
colorbar( tricontourf( Th, tau ) )
Th.plot_triangles( ax = ax, color = 'k', lw = .5, alpha = .2 )
Th.plot_boundaries( clip_on = False, color = 'grey'  )
# ax.legend()
ax.axis('equal'); ax.axis('off')
ax.set_xticks([]); ax.set_yticks([])

show()

# fespace Vh( Th, P1 );
# Vh u,v;
# ''')
#
# script += pyff.OutputScript( Th = 'mesh' )
#
# script += pyff.VarfScript(
#     stiffness = 'int2d(Th)( dx(u)*dx(v) +  dy(u)*dy(v) )',
#     Grammian = 'int2d(Th)( u*v )',
#     boundary_Grammian = 'int1d(Th, 1, 2)( u*v )'
#     )
#
# ff_output = script.get_output()
