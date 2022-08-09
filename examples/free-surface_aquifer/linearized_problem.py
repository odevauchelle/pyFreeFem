from pylab import *
from copy import deepcopy
from scipy.sparse.linalg import spsolve
from matplotlib.tri import LinearTriInterpolator
import hashlib

import style as st

import sys
sys.path.append('/home/olivier/git/pyFreeFem')

import pyFreeFem as pyff

R = 0.15
H = .8

#########################
#
# figure
#
#########################

fig, ax = subplots( ncols = 2, figsize = (8,5) )

fig = dict( z = fig, omega = fig )
ax = dict( zip( list( fig.keys() ), ax ) )

#########################
#
# build mesh
#
#########################

seepage_bottom = ( 0, 0 )
seepage_top = ( R, 0 )
divide_top = ( 1, 0 )
divide_bottom = ( 1, -H  )
river_bottom = ( 0, divide_bottom[1] )

boundary_names = [ 'free_surface', 'seepage_face',  'river_wall', 'bottom', 'divide' ]

boundaries = dict( zip( boundary_names, [
    [ divide_top, seepage_top ],
    [ seepage_top, seepage_bottom ],
    [ seepage_bottom, river_bottom ],
    [ river_bottom, divide_bottom ],
    [ divide_bottom, divide_top ],
    ] ) )

all_points = []

for boundary_name in boundary_names :

    boundary = boundaries[ boundary_name ]

    all_points += boundary[:-1]

Th = pyff.TriMesh( *array(all_points).T )

point_indices = [ 0, 0 ]

for boundary_name in boundary_names :

    point_indices[1] += len( boundaries[ boundary_name ] ) - 1

    indices = arange( point_indices[0], point_indices[1] + 1 )

    if indices[-1] >= len( all_points ) :
        indices[-1] = 0

    print( boundary_name, point_indices, indices )

    Th.add_boundary_edges( indices, boundary_name )

    point_indices[0] = point_indices[1]

error_field = Th.x*0 + 1

for _ in range(6) :

    Th = pyff.adaptmesh( Th, error_field, hmax = H/15, iso = 1, err = 1e-2 )

    #########################
    #
    # FE matrices
    #
    #########################

    script = pyff.InputScript( Th = Th )

    script += pyff.edpScript('fespace Vh( Th, P1 );')

    FE_matrices = dict( stiffness = 'int2d(Th)( dx(u)*dx(v) +  dy(u)*dy(v) )' )

    boundary_name_to_int = Th.get_boundary_label_conversion()[0]
    # print( Th.get_boundary_edges() )

    for boundary_name in boundary_names :
        FE_matrices[ boundary_name ] = 'int1d(Th,' + str( boundary_name_to_int[ boundary_name ] ) + ')( u*v )'


    script += pyff.VarfScript( **FE_matrices )
        # BoundaryGramian = 'int1d(Th, 1, 2, 3, 4 )( u*v )'

    FE_matrices = script.get_output()

    ##########################
    #
    # Find the origin
    #
    ##########################

    i_origin = set( Th.get_boundaries()['seepage_face'][0] ) & set( Th.get_boundaries()['river_wall'][0] )

    if len( i_origin ) == 1 :
        i_origin = list( i_origin )[0]
    else :
        print("Can't find origin.", i_origin )
        break

    #########################
    #
    # Solve problems
    #
    #########################

    epsilon = 1e-6

    #########################
    #
    # x
    #
    #########################

    M = - FE_matrices['stiffness']
    B = Th.x*0

    boundary_name = 'seepage_face'
    M += 1/epsilon*FE_matrices[boundary_name]

    boundary_name = 'river_wall'
    M += 1/epsilon*FE_matrices[boundary_name]

    boundary_name = 'divide'
    M += 1/epsilon*FE_matrices[boundary_name]
    B += 1/epsilon*FE_matrices[boundary_name]*( Th.x*0 + 1 )

    boundary_name = 'bottom' # x = u
    M += 1/epsilon*FE_matrices[boundary_name]
    B += 1/epsilon*FE_matrices[boundary_name]*( Th.x )

    boundary_name = 'free_surface'
    M += 1/epsilon*FE_matrices[boundary_name]
    B += 1/epsilon*FE_matrices[boundary_name]*( Th.x - R )/( 1 - R )

    x = spsolve( M, B )

    #########################
    #
    # y
    #
    #########################

    M = - FE_matrices['stiffness']
    B = Th.x*0

    boundary_name = 'bottom' # dy/dn = 1
    B += FE_matrices[boundary_name]*( Th.x*0 + 1 )

    boundary_name = 'free_surface'
    B += -FE_matrices[boundary_name]*( Th.x*0 + 1 )/( 1 - R )

    y = spsolve( M, B )

    #########################
    #
    # Error
    #
    #########################

    error_field = ( x - Th.x )*( y - Th.y )

##########################
#
# Mapping
#
##########################

Th = dict( omega = Th, z = deepcopy(Th) )

Th['z'].x = x
Th['z'].y = y
# Th['z'].y -= Th['z'].y[i_origin]

omega = Th['omega'].x + 1j*Th['omega'].y
z = Th['z'].x + 1j*Th['z'].y

Phi = 1j*( omega - z )

##########################
#
# Plots & show
#
##########################

titles = dict( z = '$z$', omega = r'$\omega$' )

for key in ax.keys() :

    ax[key].tricontour( Th[key], real(Phi), **st.flow['iso_head'] )
    ax[key].tricontour( Th[key], imag(Phi), **st.flow['flow_lines'] )

    # Th[key].plot_triangles( ax = ax[key], **st.mesh )

    for name, segments in Th[key].get_boundaries().items() :

        Th[key].plot_boundaries( { name : segments },  ax = ax[key], clip_on = False, **st.boundary[name] )

    ax[key].set_title( titles[key] )

    ax[key].axis('scaled')
    ax[key].axis('off')

# ax['omega'].legend( loc = 'center' )

show( block = True )
