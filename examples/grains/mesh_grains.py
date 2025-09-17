from pylab import *
import triangle as tr
import json as json
from scipy.sparse.linalg import spsolve
# from shapely.geometry import Polygon as shp_Polygon
from shapely import wkt as shp_wkt


with open('grains.json') as the_file :
    boundaries = json.load(the_file)

Ts = { 'segments' : [], 'vertices' : [], 'segment_markers' : [], 'holes' :[] }


for interior in shp_wkt.loads( boundaries['box'] ).interiors : # there is only one interior

    for point in interior.coords[:-1] :

        Ts['vertices'] += [ point ]

        try :
            vertex_index += 1
        except :
            vertex_index = 0

        if vertex_index > 0 :
            Ts['segments'] += [ [ vertex_index - 1, vertex_index ] ]

    Ts['segments'] += [ [ vertex_index, 0 ] ]

Ts['segment_markers'] = ( -arange( 1, 5 ) ).tolist()

shrink_factor = .97

for grain in boundaries['grains'] :

    points = []

    grain = shp_wkt.loads( grain )
    center = grain.centroid
    Ts['holes'] += [ [center.x, center.y] ]

    for point in grain.exterior.coords[:-1] :

        point = center.x + shrink_factor*( point[0] - center.x ), center.y + shrink_factor*( point[1] - center.y )

        points += [ point ]
        vertex_index += 1

        if len(points) > 1 :
            Ts['segments'] += [ [ vertex_index - 1, vertex_index ] ]
            Ts['segment_markers'] += [ -5 ]

    Ts['segments'] += [ [ vertex_index, vertex_index - len(points) + 1 ] ]
    Ts['segment_markers'] += [ -5 ]

    Ts['vertices'] += points

###################
#
# Check
#
###################

fig_grains = figure()
ax_grains = gca()
# print(len(Ts['vertices']))


for vertex in Ts['vertices'] :
    ax_grains.plot(*vertex, 'ko')

for hole in Ts['holes'] :
    ax_grains.plot(*hole, '+m')

for segment in Ts['segments'] :
    ax_grains.plot( *array(Ts['vertices'])[segment].T, color = 'm'  )

axis('scaled')

###################
#
# Build triangulation
#
###################


T = tr.triangulate( Ts, 'pa' )
tr.compare( plt, Ts, T )
fig_tri = gcf()


####################
#
# pyFreefem mesh
#
####################


import sys
sys.path.append('./../../')

import pyFreeFem as pyff

Th = pyff.triangle_to_TriMesh( T )
Th.rename_boundary( { -1 :'bottom', -2:'side', -3:'top', -4:'side', -5:'grain' } )

boundary_number, boundary_name = Th.get_boundary_label_conversion()

Th = pyff.adaptmesh( Th, hmax = .1 )

####################
#
# pyFreefem matrix
#
####################

script = pyff.InputScript( Th = 'mesh' )
script += '''
fespace Vh( Th, P2 );
Vh u,v;
Vh X = x, Y = y;

fespace VhDiff( Th, P0 );
'''
script += pyff.VarfScript(
    stiffness = 'int2d(Th)( dx(u)*dx(v) +  dy(u)*dy(v) )',
    top = 'int1d(Th,' + str( boundary_number['top'] ) + ')( u*v )',
    bottom = 'int1d(Th,' + str( boundary_number['bottom'] ) + ')( u*v )',
    dx = 'int2d(Th)( dx(u)*v )',
    dy = 'int2d(Th)( dy(u)*v )',
    )

script += pyff.VarfScript(
    gramian = 'int2d(Th)( du*v )',
    functions = ('du','v'),
    fespaces = ('VhDiff', 'Vh')
    )

script += pyff.OutputScript( X = 'vector', Y = 'vector' )

epsilon = 1e-6

for _ in range(3) :

    try :
        Th = pyff.adaptmesh( Th, c, hmin = 0.02 )
    except :
        pass

    FE_matrices = script.get_output( Th = Th )
    M = FE_matrices['stiffness'] - 1/epsilon*( FE_matrices['top'] + FE_matrices['bottom'] )
    B = FE_matrices['top']@( FE_matrices['X']*0 + 1 )
    c = spsolve( M, B )

p = pyff.get_projector( Th, 'P2', 'P1' )
gradient = spsolve( FE_matrices['gramian'].T@FE_matrices['gramian'], FE_matrices['gramian'].T@FE_matrices['dx']@c ), spsolve(FE_matrices['gramian'].T@ FE_matrices['gramian'], FE_matrices['gramian'].T@FE_matrices['dy']@c )

Xt = [] # centroids of triangles
At = [] # areas of triangles

for tr in Th.triangles :
    X_local = array( [ Th.x[tr], Th.y[tr] ] ).T
    Xt += [ mean(X_local, axis = 0) ]
    At += [ cross( X_local[1] - X_local[0], X_local[2] - X_local[0] )/2 ]

# figure()
# ax_area = gca()
# ax_area.hist( array(At)*len(Th.triangles) )
# print( 'Total area:', sum(At) )

fig_c = figure()
ax_c = gca()
ax_c.tricontourf( Th, p@c )
# ax_c.quiver( *array(Xt).T, *gradient, color = 'w'  )
contours = ax_c.tricontour( Th, p@c, colors = 'w', linestyles = '-', linewidths = .75, levels = 50 )

fig_h = figure()
ax_h = gca()
selection = abs( array( Xt )[:,1] ) < .15
a = angle( 1j*gradient[0] - gradient[1] )[ selection ]*180/pi
weights = array(At)[ selection ]
ax_h.hist( a, weights = weights, alpha = .5, bins = 150*linspace(-1,1,30), density = True )
ax_h.set_xlabel('Gradient orientation [deg]')
ax_h.set_ylabel('Probability density')

x_lim = ax_h.get_xlim()
x_fit = linspace( *x_lim, 101 )
sigma = sqrt( sum( a**2*weights )/sum(weights) )

# pdf_fit = 1/sigma*exp( - 2*abs(x_fit)/sigma )
pdf_fit = 1/(sigma*sqrt(2*pi))*exp( - (x_fit/sigma)**2/2 )
ax_h.plot(x_fit, pdf_fit, '--')
print('std angle:', sigma,  'deg')

# Th.plot_triangles( color = 'w', alpha = .1, lw = .75)#( labels = 'label' )
# Th.plot_boundaries( label = 'grains', color = 'grey' )
# Th.plot_nodes( labels = 'label', color = 'tab:blue' )
# legend()

ax_c.axis('equal')
ax_c.axis('off')

# fig_h.savefig( 'grains_histogram.svg', bbox_inches = 'tight' )
# fig_c.savefig( 'grains_field.svg', bbox_inches = 'tight' )
# fig_tri.savefig( 'grains_compare.svg', bbox_inches = 'tight' )
# fig_grains.savefig( 'grains_mesh.svg', bbox_inches = 'tight' )


show()
