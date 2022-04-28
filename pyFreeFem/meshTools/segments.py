
import numpy as np
import warnings
from itertools import permutations

def concatenate_segments( tail, head) :
    '''
    body = concatenate_segments( tail, head )

    Parameters:
        tail (list) : list of integers
        head (list) : list of integers

    Returns:
        body(list) or None: concatenated head and tail if they match
    '''

    if tail[-1] == head[0]:
        return tail + head[1:]


def edges_to_segments( edges ) :
    '''
    segments = edges_to_segments( edges )

    Assemble individual edges into boundary segments. Might be unstable if nodes (integer) are repeated
    in edges more than twice, that is, if the boundary is branching.

    Parameters :
        edges (list) : list of edges. An edge is a list of two nodes (integers).

    Returns :
        segments (list) : list of segments. A segment is a list of nodes.
    '''

    segments = edges.copy()

    keep_going = True

    while keep_going :

        keep_going = False

        for perm in permutations( range( len( segments ) ), 2 ) :

            new_segment = concatenate_segments( segments[perm[0]], segments[perm[1]] )

            if not new_segment is None :

                for i in sorted( perm, reverse = True ) :
                    segments.pop( i )

                segments += [ new_segment ]

                keep_going = True

                break

    return segments

def label_conversion( labels ) :

    kept_labels = []
    changed_labels = []

    for label in list( set( labels ) ) :

        if isinstance( label, ( int, np.int64 ) ) :
            kept_labels += [ label ]

        else :
            changed_labels += [ label ]

    # make sure order is consistent; not guaranteed by set
    changed_labels.sort( key = lambda label: ( type( label ).__name__, label ) )

    new_int_label_start = invent_label( kept_labels )

    new_int_labels = list( range( new_int_label_start, new_int_label_start + len( changed_labels ) ) )

    int_labels = kept_labels + new_int_labels
    labels = kept_labels + changed_labels

    return dict( zip( labels, int_labels ) ), dict( zip( int_labels, labels ) )


def triangle_edge_to_node_edge( triangle_edge, triangles ) :
    '''
    ( triangle_index, node_index_in_triangle ) -> ( start_node_index, end_node_index )
    '''

    triangle_index, triangle_node_index = triangle_edge

    start_node_index = triangles[ triangle_index, triangle_node_index ]

    if triangle_node_index < 2 :
        triangle_end_node_index = triangle_node_index + 1
    else :
        triangle_end_node_index = 0

    end_node_index = triangles[ triangle_index, triangle_end_node_index ]

    return start_node_index, end_node_index

# triangle_edge_to_node_edge = np.vectorize(triangle_edge_to_node_edge)

def find_triangle_index( triangles, start_node, end_node ) :
    '''
    find_triangle_index( triangles, start_node, end_node )
        triangles : triangulation
        start_node, end_node : nodes indices of edge

    returns the triangle_index of this edge. There can only be one, since orientation is strict.
    '''

    triangle_index = None

    for nodes_order in [ [0,1], [1,2], [-1,0] ] :
        try :
            triangle_index = triangles[ :, nodes_order].tolist().index( [ start_node, end_node ] )
            break
        except :
            pass

    return triangle_index

def edge_nodes_to_triangle_edge( edge_nodes, triangles, flip_reversed_edges = True ) :
    '''
    ( start_node, end_node ) -> ( triangle_index, start_node_in_triangle )
    '''

    start_node, end_node = edge_nodes

    triangle_index = find_triangle_index( triangles, start_node, end_node )

    if triangle_index is None and flip_reversed_edges:
        warnings.warn('Reversing some edges')
        start_node, end_node = end_node, start_node
        triangle_index = find_triangle_index( triangles, start_node, end_node )

    if triangle_index is None :
        warnings.warn('Could not find some boundary edges. They are lost.' )
        return None

    else :
        start_node_in_triangle = triangles[triangle_index].tolist().index( start_node )
        return triangle_index, start_node_in_triangle

def nodes_to_edges( nodes, label = None ) :
    '''
    node_indices, label -> [ [ start_node_index, end_node_index, label ], ... ]
    '''
    return list( zip( *[ nodes[:-1], nodes[1:], [label]*( len( nodes ) - 1 ) ] ) )

def edges_to_boundary_edges( edges ) :
    '''
    start_node, end_node, label -> { (start_node, end_node) : label }
    or
    triangle, start_node, label -> { (triangle, start_node) : label }
    '''

    edges_dict = {}

    for edge in edges :
        edges_dict.update( { tuple( edge[:-1] ) : edge[-1] } )

    return edges_dict

def invent_label( existing_labels ) :
    existing_int_labels = [ label for label in existing_labels if type(label) is int ]
    return max( existing_int_labels + [0] ) + 1


def node_index_to_triangle_index_edges( edges, triangles, label = 1 ) :
    '''
    triangle_index_edges = edges_node_index_egdes_to_triangle_edges( edges, triangles, label = 1 )

    Convert a list of segments indexed by node to a list of segment indexed by triangle.

    Arguments:
    edges : a list of segments indexed by node ( [ [i, j], ... ], where i and i are node indices )
    triangles : a triangulation
    label : the label of the edges
    '''
    return [ list( edge_nodes_to_triangle_edge( edge, triangles ) ) + [ label ] for edge in edges ]


if __name__ == '__main__' :

    from random import shuffle

    # from pylab import *
    # import matplotlib.tri as tri
    #
    # theta =  linspace( 0, 2*pi, 6 )
    # x, y = cos(theta), sin(theta)
    # x, y = append(x,0), append(y,0)
    #
    # triangles = array( [ [i, i + 1, len(x)-1] for i in range(len(x)-1) ] )
    #
    # edge_nodes = range(len(x)-1)
    # label = 'toto'
    #
    # print(triangles)
    #
    # mesh = tri.Triangulation( x, y, triangles = triangles )
    #
    # triplot(mesh)
    # plot( x[edge_nodes], y[edge_nodes] )
    #
    # print( edges_to_boundary_edges( nodes_to_edges( edge_nodes, label ) ) )
    # labels = ['toto','tortue', 1, 6, '3', 45.6, 12, 'chat', 35, '54' ]
    labels = [5,'5',12,12.3]

    print( invent_label( labels ) )

    # for table in label_conversion( labels )[0] :
    print('------')

    tables = []

    for _ in range(10) :
        tables += [ label_conversion( labels )[1] ]
        shuffle( labels )

    for key in tables[0].keys() :

        print( (key,), ':', *[ (table[key],) for table in tables ] )


    # axis('equal')
    # show()
