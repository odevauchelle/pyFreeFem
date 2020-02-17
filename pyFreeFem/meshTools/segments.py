######################

import numpy as np

def find_disjunctions( seg_list ) :
    '''
    Find disjunctions (holes) in a list of segments.

    Parameters :
    ----------
    seg_list : numpy.array of segments
        list of couples of nodes coordinates.

    Returns :
    ----------
    disjuctions : numpy.array of indices
        Indices of the disjunction segments

    Example :
    ----------
    disjuctions = find_disjunction( array( [ [1,12], [12,3], [3,7], [4, 14], [14,1] ] ) )
    >>> [3]

    '''

    if len( seg_list ) < 2 :
        return seg_list

    else :
        return ( seg_list[ 1:-1, 0 ] - seg_list[ 0:-2, 1 ] ).nonzero()[0] + 1

def stitch_segment_list( seg_list, n, verbose = False ) :

    '''
    Cut a segment list at location 'n', finds the next occurence of the end node, and stich the rest of the list there.

    Parameters :
    ----------
    seg_list : numpy.array of segments
        List of couples of nodes coordinates.
    n : integer
        Location of the disjuction to be mended.

    Returns :
    ----------
    mended_segment_list : numpy.array of indices
        Indices of the disjunction segments
    mending : Boolean
        Whether the list of segment was mended

    Example :
    ----------
    stitch_segment_list( array( [ [1,12], [12,3], [7, 14], [3, 7], [14, 1 ] ] ), 2 )
    >>> (array([[ 1, 12], [12,  3], [ 3,  7], [14,  1], [ 7, 14]]), True)
    '''

    try :
        n2 = n + np.where( seg_list [ n:, 0  ] == seg_list[ ( n - 1 ), 1 ] )[0][0]
        return np.concatenate( ( seg_list[ 0 : n ], seg_list[ n2 : ], seg_list[ n : n2 ] )  ), True

    except :
        if verbose :
            print('No possible reconnection of segment list. Returning initial segment list.')
        return seg_list, False

def reorder_boundary( seg_list ) :
    '''
    Reorder a list of segments until there is no disjunction left, or until the remaining disjunctions cannot be solved.

    Parameters :
    ----------
    seg_list : numpy.array of segments
        List of couples of nodes coordinates.

    Returns :
    ----------
    ordered_segment_list : numpy.array of indices
        List of ordered segments, with the least possible disjuctions

    Example :
    ----------
    reorder_boundary( [ [1,12], [12,3], [7, 14], [3, 7], [14, 1 ] ] )
    >>> [[ 1 12], [12  3], [ 3  7], [ 7 14], [14  1]]
    '''

    keep_trying = True

    seg_list = np.array(seg_list)

    while keep_trying :

        keep_trying = False

        for disjunction_index in find_disjunctions( seg_list ) :

            seg_list, mending = stitch_segment_list( seg_list, disjunction_index )

            if mending :
                keep_trying = True
                break

    return seg_list.tolist()

def ordered_edges_to_segments( edges ) :

    segments = []
    segment = edges[0]
    previous_edge = edges[0]

    for edge in edges[1:] :

        if edge[0] == previous_edge[1] :
            segment += [ edge[1] ]

        elif len(segment) > 0 :
            segments += [ segment ]
            segment = edge

        previous_edge = edge

    if len(segment) > 0 :
        segments += [ segment ]

    return segments

def label_conversion( labels ) :

    labels = list( set( labels ) )
    int_labels = range( 1, len( labels ) + 1 )

    int_to_label, label_to_int =  dict( zip( labels, int_labels ) ), dict( zip( int_labels, labels ) )

    return int_to_label, label_to_int


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

if __name__ == '__main__' :

    labels = ['toto', 2, 1]

    none_and_int = label_conversion( ['toto', 2, 1] )
    print( none_and_int )

    ordered_edges = reorder_boundary( [ [1,12], [12,3], [7, 14], [8, 7], [14, 20], [7, 9], [12, 9] ] )
    print( ordered_edges  )
    print( ordered_edges_to_segments( ordered_edges  ) )
