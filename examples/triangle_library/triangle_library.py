from pylab import *
import triangle as tr

'''
https://rufat.be/triangle/index.html
'''

Ts = dict(
    vertices = [ [-1,-1], [1,-1], [1,1], [-1,1], [0.6,-1], [-0.8, 0], [-0,-.1] ],
    segments = [ [ 0, 4 ], [ 4, 1 ], [1, 2], [2,3], [3, 0], [4,5], [5, 6] ],
    )

Ts['segment_markers'] =  [-1]*5 + [-2]*(len(Ts['segments']) - 5 )

T = tr.triangulate( Ts, 'pa' )

for key in T.keys() :
    print(key)

tr.compare(plt, Ts, T)
savefig( 'compare.svg', bbox_inches = 'tight' )



import sys
sys.path.append('./../../')

import pyFreeFem as pyff

Th = pyff.triangle_to_TriMesh( T )

figure()

Th.rename_boundary( { -1 :'box', -2:'river' } )

for _ in range(3):
    Th = pyff.adaptmesh( Th, hmax = .2, iso = 1 )

Th.plot_triangles()
Th.plot_boundaries()
legend()

axis('equal')
axis('off')
savefig( 'river.svg', bbox_inches = 'tight' )

show()
