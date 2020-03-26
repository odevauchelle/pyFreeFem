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



edges = array( list(  Th.boundary_edges.keys() ) )[ array( list(  Th.boundary_edges.values() ) ) == 3 ]
edges = [ pyff.triangle_edge_to_node_edge( e, Th.triangles ) for e in edges  ]
edges = [ [1,12], [12,3], [7, 14], [3, 7], [14, 1 ] ]
# edges = edges[::-1]

def reorder_boundary( edges ):
    '''
    ordered_edges = reorder_boundary( edges )

    Arguments:
        edges (list) : list of edges. Each edge is a tuple of two node indices.
    '''

    segments = []

    while len(edges) > 0 :

        edge = edges.pop(0)

        for segment in segments :

            if segment[-1][-1] == edge[0] :
                segment += [ edge ]
                edge = []
                break

            elif segment[0][0] == edge[-1] :
                segment = [ edge ] + segment
                edge = []
                break

        if edge != [] :
            segments += [ [ edge ] ]

    return segments



print('edges')
print(edges)
print(pyff.reorder_boundary(edges))
print(reorder_boundary(edges))

#
# print('Th.get_boundary_edges()')
# print( Th.get_boundary_edges() )
#
#
#
# #
# # plot( Th.x[boundary_nodes], Th.y[boundary_nodes], color = 'red' )
# #
#
# ax.legend()
# ax.axis('equal'); ax.axis('off')
# ax.set_xticks([]); ax.set_yticks([])
#
# show()
