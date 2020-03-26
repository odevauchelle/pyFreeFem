from pylab import *
from scipy.sparse.linalg import spsolve

import sys
sys.path.append('./../../')

import pyFreeFem as pyff

script = pyff.edpScript('mesh Th = square( 3, 3);')
script += pyff.OutputScript( Th = 'mesh' )
Th = script.get_output()['Th']


figure()
ax = gca()

#colorbar( tricontourf( Th, tau ) )
Th.plot_boundaries( clip_on = False )
Th.plot_triangles( ax = ax, labels = 'index', color = 'k', lw = .5, alpha = .2 )
Th.plot_nodes( labels = 'index', color = 'grey' )

print('Th.boundary_edges')
print(Th.boundary_edges)

print('Th.get_boundary_edges()')
print( Th.get_boundary_edges() )

#
# plot( Th.x[boundary_nodes], Th.y[boundary_nodes], color = 'red' )
#

ax.legend()
ax.axis('equal'); ax.axis('off')
ax.set_xticks([]); ax.set_yticks([])

show()
