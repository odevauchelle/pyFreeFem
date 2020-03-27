


edges = [ [ 6, 7 ], [ 5, 6 ], [ 7, 12 ], [ 15, 18 ], [ 12, 15 ], [34,67], [876,34], [ 7, 14 ] ]

from numpy.random import shuffle

for _ in range(10):
    shuffle(edges)
    print(*edges)
    print(*edges_to_segments(edges))
    print('------------------------')
