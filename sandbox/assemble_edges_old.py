from random import shuffle

def concatenate_segments( segment_1, segment_2 ) :

    if segment_1[-1][1] == segment_2[0][0] :
        return [ segment_1 + segment_2 ]

    elif segment_2[-1][1] == segment_1[0][0] :
        return [ segment_2 + segment_1 ]

    else :
        return [ segment_2, segment_1 ]

edges = [ [ 6, 7 ], [ 5, 6 ], [ 7, 12 ], [ 15, 18 ], [ 12, 15 ], [ 6, 7 ] ]

segments = [ [ edge ] for edge in edges ]


print(segments)

for _ in range(5) :

    segments = concatenate_segments( segments.pop(0), segments.pop(0) ) + segments
    shuffle(segments)
    print(len(segments))
