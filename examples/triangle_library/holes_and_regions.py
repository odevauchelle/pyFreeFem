from pylab import *
import triangle as tr

sys.path.append('./../../')
import pyFreeFem as pyff


theta = linspace(0,2*pi,12)[1:]

vertices = array( [ cos(theta), sin(theta) ] ).T
segments = list( array( [ arange( len( vertices ) - 1 ), arange( 1, len(vertices) ) ] ).T ) + [ [ len(vertices)-1, 0 ] ]

vertices = list( vertices ) + list( vertices*.5 )
segments = segments + list( array( segments ) + len( segments ) )

T_hole = tr.triangulate( dict( vertices = vertices, segments = segments, holes = [ [ 0, 0 ] ] ), 'pa' )
T =  tr.triangulate( dict( vertices = vertices, segments = segments, regions = [ [ 0, 0, -1, 0 ] ] ), 'pA' )

for key in T.keys():
    print(key)

print( T['segment_markers'] )

tr.compare( plt, T_hole,T ) 
savefig('hole.svg',bbox_inches = 'tight')


Th = pyff.triangle_to_TriMesh( T )

figure()
Th.plot_triangles( labels = 'label' )
Th.plot_boundaries()
Th.plot_nodes( labels = 'label', color = 'tab:blue' )
legend()

print(Th.boundary_edges)

axis('off')
axis('equal')

savefig('region.svg',bbox_inches = 'tight')


show()