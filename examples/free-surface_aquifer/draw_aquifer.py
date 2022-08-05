from pylab import *
import style as st

H = .3

boundaries = dict(
    bottom = -1j*H + array( [ 0, 1 ] ),
    river_wall = 1j*array( [ -H, 0] ),
    free_surface = array( [0] + list( logspace( -3,0,30 ) ) ) + 0*1j
)

boundaries['free_surface'] += 0.2*1j*( sqrt( 1 - ( boundaries['free_surface'] - 1 )**2 ) + 1 )

boundaries['seepage_face'] = array([0,1])*boundaries['free_surface'][0]

boundaries['divide'] = array([boundaries['bottom'][-1], boundaries['free_surface'][-1]])

######################
#
# Plot boundaries
#
######################

ax = gca()

st.plot_boundaries( boundaries, ax )

ax.set_xlabel('$x$')
ax.set_ylabel('$y$')

ax.text( .5, mean( imag(boundaries['bottom'])), '$z=x+iy$\n', ha = 'center', va = 'bottom' )

####################
#
# Arrows
#
####################

n_rain_arrow = 5
dx_rain = 1/n_rain_arrow
x_rain = linspace( -dx_rain/2, 1 + dx_rain/2, n_rain_arrow + 2 )[1:-1]
y_rain = imag( boundaries['free_surface'][-1] ) + 0.15

for x in x_rain :
    ax.add_patch( st.arrow( x, y_rain ) )

ax.add_patch( st.arrow( 0.0, imag(mean(boundaries['seepage_face'])), orientation = 'right' ) )

ax.set_ylim( ax.get_ylim()[0], y_rain )

####################
#
# Filling
#
####################

st.fill_in( boundaries, ax )

# savefig('../../figures/aquifer_boundaries.svg', bbox_inches = 'tight')

show()
