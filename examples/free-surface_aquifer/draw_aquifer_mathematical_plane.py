from pylab import *
import style as st

Q = .2

boundaries = dict(
    seepage_face = array([0,Q]) + 0j,
    free_surface = array([Q,1]) + 0j,
    bottom = linspace( 0, 1, 30 ) + 0j,
)


boundaries['bottom'] += cos( pi*boundaries['bottom'] )*.05j - .6j

boundaries['divide'] = array([boundaries['bottom'][-1], boundaries['free_surface'][-1]])
boundaries['river_wall'] = array([boundaries['bottom'][0], boundaries['seepage_face'][0]])

######################
#
# Plot boundaries
#
######################

ax = gca()

st.plot_boundaries( boundaries, ax )

ax.set_xlabel('$u$')
ax.set_ylabel('$v$')

ax.text( .5, mean( imag(boundaries['bottom'])), '$\omega=u+iv$\n', ha = 'center', va = 'bottom' )

ax.plot( Q, 0, 'ok', color = st.boundary['free_surface']['color'])
ax.text( Q, 0, '\n$\omega=R$', ha = 'center', va = 'top', color = st.boundary['free_surface']['color'])

####################
#
# Filling
#
####################

st.fill_in( boundaries, ax )

savefig('../../figures/aquifer_boundaries_mathematical_plane.svg', bbox_inches = 'tight')

show()
