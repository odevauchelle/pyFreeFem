from pylab import *
import triangle as tr
import json as json
# from shapely.geometry import Polygon as shp_Polygon
from shapely import wkt as shp_wkt

'''
https://rufat.be/triangle/index.html
'''

with open('grains.json') as the_file :
    boundaries = json.load(the_file)


Ts = { 'segments' : [], 'vertices' : [], 'segment_markers' : [] }


for interior in shp_wkt.loads( boundaries['box'] ).interiors : # there is only one interior

    for point in interior.coords :

        Ts['vertices'] += [point]

        try :
            vertex_index += 1
        except :
            vertex_index = 0

        if vertex_index > 0 :
            Ts['segments'] += [ [ vertex_index - 1, vertex_index ] ]
            Ts['segment_markers'] += [ -1 ]



for grain in boundaries['grains'] :

    points = []

    for point in shp_wkt.loads( grain ).exterior :
        points += [ point ]
        vertex_index += 1
        Ts['segments'] += [ [ vertex_index - 1, vertex_index ] ]
        Ts['segment_markers'] += [ -2 ]

    Ts['segments'] += [ [ vertex_index, vertex_index - len(points) ] ]
    Ts['segment_markers'] += [ -2 ]

    Ts['vertices']



# print(Ts)

#
# Ts = dict(
#     vertices = [ [-1,-1], [1,-1], [1,1], [-1,1], [0.6,-1], [-0.8, 0], [-0,-.1] ],
#     segments = [ [ 0, 4 ], [ 4, 1 ], [1, 2], [2,3], [3, 0], [4,5], [5, 6] ],
#     )
#
# Ts['segment_markers'] =  [-1]*5 + [-2]*(len(Ts['segments']) - 5 )

T = tr.triangulate( Ts )

for key in T.keys() :
    print(key)

tr.compare(plt, Ts, T)

show()
# savefig( 'compare.svg', bbox_inches = 'tight' )

#
# import sys
# sys.path.append('./../../')
#
# import pyFreeFem as pyff
#
# Th = pyff.triangle_to_TriMesh( T )
#
# figure()
#
# Th.rename_boundary( { -1 :'box', -2:'river' } )
#
# for _ in range(3):
#     Th = pyff.adaptmesh( Th, hmax = .2, iso = 1 )
#
# Th.plot_triangles()
# Th.plot_boundaries()
# legend()
#
# axis('equal')
# axis('off')
# # savefig( 'river.svg', bbox_inches = 'tight' )

# show()
