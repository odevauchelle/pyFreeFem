
from pylab import *
import matplotlib.tri as tri

################# PARAMETERS
L = 1.
H = .5
k = pi/L
levels = logspace( -3, 0, 15 )

def Phi(z) :
    return -cosh( k*( 1j*z - H ) )

################ MESH
x = linspace( 0, L/2., 50 )
y = linspace( -H, 0, int( len(x)/max(x)*H ) )
x, y = meshgrid( x, y )
x = x.flatten(); y = y.flatten()
z = x + 1j*y

Th = tri.Triangulation( x, y )

################# PLOT
triplot(Th, lw = .7, color = 'k', alpha = .2)
contours = tricontour( Th, imag( Phi(z) ), colors = ['tab:blue'], levels = levels ).collections
axis('equal'); axis('off')

################# TRAVEL TIMES

x_start = []
t = []

for contour in contours:
    for path in contour.get_paths() :
        x, y = path.vertices.T
        z = x + 1j*y
        phi = real( Phi(z) )
        ds = sqrt( diff(x)**2 + diff(y)**2 )
        dt = ds**2/diff(phi)
        x_start += [ x[0] ]
        t += [ sum(dt) ]

figure()
plot( x_start, t )


show()
