from pylab import *


H = .3


boundaries = dict(
    bottom = -1j*H + array( [ 0, 1 ] ),
    river_wall = 1j*array( [ -H, 0] ),
    free_surface = linspace(0,1,50) + 0*1j
)

boundaries['free_surface'] += 0.2*1j*( sqrt( 1 - ( boundaries['free_surface'] - 1 )**2 ) + 1 )

boundaries['seepage face'] = array([0,1])*boundaries['free_surface'][0]

boundaries['divide'] = array([boundaries['bottom'][-1], boundaries['free_surface'][-1]])

for name, z in boundaries.items() :
    color = plot( real(z), imag(z), label = name.replace('_', ' ') )[0].get_color()

legend()
axis('equal')

savefig('../../figures/aquifer_boundaries.svg', bbox_inches = 'tight')

show()
