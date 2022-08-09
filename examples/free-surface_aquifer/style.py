from matplotlib.patches import FancyArrowPatch, ArrowStyle
from pylab import *
from matplotlib.patches import Polygon

##################
#
# Arrow patches
#
##################

arrow_length = .1
arrow_style = ArrowStyle( "simple", head_length=.3 )


def arrow( x, y, orientation = 'down' ):

    X_start = array( [ x, y ] )
    X_end = X_start.copy()

    if orientation == 'down' :
        X_end[1] -= arrow_length

    elif orientation == 'right': # axis are reversed !
        X_end[0] -= arrow_length


    return FancyArrowPatch( X_start, X_end, mutation_scale = 25, edgecolor = 'none', alpha = .3, arrowstyle = arrow_style, clip_on = False )

flow = dict(
    iso_head = dict( colors = 'tab:blue', linestyles = ':', alpha = .5 ),
    flow_lines = dict( colors = 'tab:blue', alpha = .5 )
    )

boundary = dict(
    free_surface = dict( color = 'tab:blue' ),
    river_wall = dict( color = 'grey', linestyle = ':'),
    seepage_face = dict( color = 'tab:orange'),
    bottom = dict( color = 'grey', linestyle = '-' ),
    divide = dict( color = 'grey', linestyle = '--')
)

mesh = dict( lw = .5, color = 'grey', alpha = .2)


def plot_boundaries( boundaries, ax = None ) :

    if ax is None :
        ax = gca()

    for name, z in boundaries.items() :
        z = flatten( array( z ) )
        ax.plot( real(z), imag(z), label = name.replace('_', ' '), **boundary[name] )

    # ax.legend()
    ax.axis('equal')




def fill_in( boundaries, ax = None  ) :

    if ax is None :
        ax = gca()

    omega = []

    for z in boundaries.values() :
        omega += list( z )

    omega = array( list( set( omega ) ) )
    omega = omega[ argsort( angle( omega - .5 ) ) ]

    ax.add_patch( Polygon( [ [real(omega), imag(omega)] for omega in omega ], edgecolor = 'none', facecolor = 'tab:blue', alpha = .1) )
