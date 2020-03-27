from pylab import *
from scipy.sparse.linalg import spsolve

import sys
sys.path.append('./../../')

import pyFreeFem as pyff

script = pyff.edpScript('mesh Th = square( 20, 20 );')
script += pyff.OutputScript( Th = 'mesh' )
Th = script.get_output()['Th']

u = Th.x*( 1 - Th.y )**2
boundary = Th.get_boundaries()[1]




ax = gca()

ax.tricontourf( Th, u )
Th.plot_boundaries( ax = ax, clip_on = False, color = 'k' )
Th.plot_triangles( ax = ax, color = 'k', lw = .5, alpha = .2 )

ax.plot( Th.x[boundary], Th.y[boundary], 'r--', clip_on = False, lw = 2 )

# ax.legend()
ax.axis('equal'); ax.axis('off')
ax.set_xticks([]); ax.set_yticks([])

show()
