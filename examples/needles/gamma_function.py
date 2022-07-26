from pylab import *
from mpmath import gammainc

def Phi_func( z, kappa, integration_const = 0 ) :
    result = exp( 2*z/kappa )*( gammainc( 1/2, a = 0, b = 2*z/kappa ) + integration_const )
    return result


x_min = -1
x = linspace( x_min, .1, 30 )
y = linspace( -1, 1, len(x) + 1)*.4

x, y = meshgrid(x,y)
z = x + 1j*y

kappa = .6

Phi = array( [ Phi_func( the_z, kappa ) for the_z in z.flatten() ], dtype = np.complex ).reshape( shape(z) )

contourf( x, y, real( Phi ) )
plot( [x_min, 0], [0]*2, color = 'tab:orange' )

axis('scaled')
show()
