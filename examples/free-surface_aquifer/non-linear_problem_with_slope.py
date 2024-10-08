from pylab import *
from copy import deepcopy
from scipy.sparse.linalg import spsolve

import style as st

import sys
sys.path.append('/home/olivier/git/pyFreeFem')

import pyFreeFem as pyff

H = .2
epsilon = 1e-6
bottom_slopes = linspace( 0, 0.4, 10 )
Rs = linspace(1, 1, len(bottom_slopes))*0.2

#########################
#
# build mesh
#
#########################

points = dict(
    seepage_bottom = ( 0, 0 ),
    seepage_top = ( Rs[0], 0 ),
    divide_top = ( 1, 0 ),
    divide_bottom = ( 1, - H ),
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



#########################
#
# FE matrices
#
#########################

script = pyff.InputScript( Th = 'mesh' )

script += pyff.edpScript('fespace Vh( Th, P1 );')

boundary_name_to_int = Th.get_boundary_label_conversion()[0]

FE_matrices = dict(
    stiffness = 'int2d(Th)( dx(u)*dx(v) +  dy(u)*dy(v) )',
    bottom_with_slope = 'int1d(Th,' + str( boundary_name_to_int[ 'bottom' ] ) + ')( u*v/sqrt( 1 + (N.x/N.y)^2 ) )'
    )


# print( Th.get_boundary_edges() )

for boundary_name in Th.get_boundaries().keys() :
    FE_matrices[ boundary_name ] = 'int1d(Th,' + str( boundary_name_to_int[ boundary_name ] ) + ')( u*v )'

script += pyff.VarfScript( **FE_matrices )

figure()
ax_curr = gca()

for i_p, bottom_slope in enumerate( bottom_slopes ) :
    
    R = Rs[i_p]

    print( 'bottom_slope', bottom_slope )
    print( 'R', R )

    for relax_i in range( 20 ) :

        print('relax_i', relax_i)
        ax_curr.cla()
        Th.plot_triangles( ax = ax_curr, linewidth = .5 )
        ax_curr.plot( 0, 0, '.', color = 'tab:orange')
        ax_curr.axis('equal')
        pause(0.01)

        try :
            Th.y += - 0.01*delta_v
            del x
        except :
            pass

        for adapt_i in range(4) :

            try :
                Th = pyff.adaptmesh( Th, x*y, iso = 1, err = 1e-3 )
                # Th = pyff.adaptmesh( Th, ( x - Th.x )*( y - Th.y ), iso = 1 )
            except :
                Th = pyff.adaptmesh( Th, 1, hmax = H/5 )

            FE_matrices = script.get_output( Th = Th )

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

            # dy/dn = 1/sqrt( 1 + (dv/du)**2 )
            B += FE_matrices['bottom_with_slope']*( Th.x*0 + 1 )

            boundary_name = 'free_surface'
            B += -FE_matrices[boundary_name]*( Th.x*0 + 1 )/( 1 - R )

            y = spsolve( M, B )


        #################
        #
        # find the origin
        #
        ##################

        y -= max( y[ Th.get_boundaries()['river_wall'][0] ] )

        ################
        #
        # movemesh
        #
        ###############

        M = - FE_matrices['stiffness']

        boundary_name = 'bottom' # x = u
        M += 1/epsilon*FE_matrices[boundary_name]
        B += 1/epsilon*FE_matrices[boundary_name]*( y + H - bottom_slope*x )

        for boundary_name in ['free_surface', 'seepage_face' ] : # delta_v = 0
            M += 1/epsilon*FE_matrices[boundary_name]

        delta_v = spsolve( M, B )



##########################
#
# Mapping
#
##########################

Th = dict( omega = Th, z = deepcopy(Th) )

Th['z'].x = x
Th['z'].y = y

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
    # ax[space].tricontour( Th[space], imag(Phi), levels = [1e-3],colors = 'red' )

    # Th[space].plot_boundaries( ax = ax[space] )
    # Th[space].plot_triangles( ax = ax[space], **st.mesh )

    for name, segments in Th[space].get_boundaries().items() :

        Th[space].plot_boundaries( { name : segments },  ax = ax[space], clip_on = False, **st.boundary[name] )

    ax[space].set_title( titles[space] )
    #
    ax[space].axis('scaled')
    ax[space].axis('off')

ax['omega'].legend( loc = 'center' )

# savefig('../../figures/free-surface_non_linear.svg', bbox_inches = 'tight')
##########################
#
# Interpolation
#
##########################

# from matplotlib.tri import LinearTriInterpolator
# measurement_point = .5, -H/10
# print( LinearTriInterpolator( Th['z'], real( Phi ) )( *measurement_point ) )


# # p['phi'] = str( LinearTriInterpolator( Th['z'], real( Phi ) )( *p['X_measurement'] ) )
# ax['z'].plot( *measurement_point, '.k' )

show( block = True )
