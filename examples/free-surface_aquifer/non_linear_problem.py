from pylab import *
from copy import deepcopy
from scipy.sparse.linalg import spsolve


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

fig_mesh = figure()
ax_mesh = gca()

#########################
#
# build mesh
#
#########################

points = dict(
    seepage_bottom = ( 0, 0 ),
    seepage_top = ( R, 0 ),
    divide_top = ( 1, 0 ),
    divide_bottom = ( 1, -H  ),
    river_bottom = ( 0, -H )
    )

boundaries = dict(
    free_surface = [ 'divide_top', 'seepage_top' ],
    seepage_face = [ 'seepage_top', 'seepage_bottom' ],
    river_wall = [ 'seepage_bottom', 'river_bottom' ],
    bottom = [ 'river_bottom', 'divide_bottom' ],
    divide = [ 'divide_bottom', 'divide_top' ]
    )

Th = pyff.TriMesh( *array( list( points.values() ) ).T )

for boundary_name, boundary_points in boundaries.items() :
    Th.add_boundary_edges( [ list( points.keys() ).index( point ) for point in boundary_points ], boundary_name )
#
Th = pyff.adaptmesh( Th, 1, iso = 1, hmax = .3 )

Th.plot_triangles( ax = ax_mesh, **st.mesh )
# Th.plot_boundaries( ax = ax_mesh )

for name in Th.get_boundaries().keys() :
    Th.plot_boundaries( [name], ax = ax_mesh, **st.boundary[name] )

ax_mesh.legend()
ax_mesh.axis('scaled')

#########################
#
# FE matrices
#
#########################

script = pyff.InputScript( Th = 'mesh' )

script += pyff.edpScript('fespace Vh( Th, P1 );')

FE_matrices = dict( stiffness = 'int2d(Th)( dx(u)*dx(v) +  dy(u)*dy(v) )' )

boundary_name_to_int = Th.get_boundary_label_conversion()[0]
# print( Th.get_boundary_edges() )

for boundary_name in Th.get_boundaries().keys() :
    FE_matrices[ boundary_name ] = 'int1d(Th,' + str( boundary_name_to_int[ boundary_name ] ) + ')( u*v )'


script += pyff.VarfScript( **FE_matrices )

for _ in range(6) :

    try :
        Th = pyff.adaptmesh( Th, ( x - Th.x )*( y - Th.y ), hmax = H/15, iso = 1, err = 1e-2 )
    except :
        pass

    FE_matrices = script.get_output( Th = Th )

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

    # error_field = ( x - Th.x )*( y - Th.y )

##########################
#
# Mapping
#
##########################

Th = dict( omega = Th, z = deepcopy(Th) )

Th['z'].x = x
Th['z'].y = y
Th['z'].y -= max( Th['z'].y[ Th['z'].get_boundaries()['river_wall'][0] ] )

omega = Th['omega'].x + 1j*Th['omega'].y
z = Th['z'].x + 1j*Th['z'].y



Phi = 1j*( omega - z )

##########################
#
# Plots & show
#
##########################

titles = dict( z = '$z$', omega = r'$\omega$' )


fig, ax = subplots( ncols = 2, figsize = (8,5) )
ax = dict( zip( list( Th.keys() ), ax[::-1] ) )

for space in Th.keys() :

    ax[space].tricontour( Th[space], real(Phi), **st.flow['iso_head'] )
    ax[space].tricontour( Th[space], imag(Phi), **st.flow['flow_lines'] )

    # Th[space].plot_boundaries( ax = ax[space] )
    # Th[space].plot_triangles( ax = ax[space], **st.mesh )

    for name, segments in Th[space].get_boundaries().items() :

        Th[space].plot_boundaries( { name : segments },  ax = ax[space], clip_on = False, **st.boundary[name] )

    ax[space].set_title( titles[space] )
    #
    ax[space].axis('scaled')
    ax[space].axis('off')

# ax['omega'].legend( loc = 'center' )

# savefig('../../figures/free-surface_non_linear.svg', bbox_inches = 'tight')

show( block = True )
