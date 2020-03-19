
from pylab import *
import matplotlib.tri as tri

################# PARAMETERS
L = 1.
H = .3
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

# fig, (ax_time, ax_psi) = subplots( nrows = 2, sharex = True, figsize = (6,6))
fig = plt.figure(figsize = (6,5))
gs = fig.add_gridspec( 3, 1 )
ax_time = fig.add_subplot( gs[2:, :] )
ax_psi = fig.add_subplot( gs[:2, :], sharex = ax_time)

ax = ax_psi
ax.triplot(Th, lw = .7, color = 'k', alpha = .2)
contours = ax.tricontour( Th, imag( Phi(z) ), colors = ['tab:blue'], levels = levels )
ax.axis('equal'); ax_psi.axis('off')

################# TRAVEL TIMESx

x_start = []
t = []

for contour in contours.collections:
    for path in contour.get_paths() :
        x, y = path.vertices.T
        z = x + 1j*y
        phi = real( Phi(z) )
        ds = sqrt( diff(x)**2 + diff(y)**2 )
        dt = ds**2/diff(phi)
        x_start += [ x[0] ]
        t += [ sum(dt) ]


trav_time_color = 'tab:red'
ax_time.plot( x_start, t, color = trav_time_color )
ax_time.set_xlabel( 'Starting position' )
ax_time.set_ylabel( 'Travel time' )

fig_path_and_name = './../../figures/' + __file__.split('/')[-1].split('.')[0] + '.svg'
savefig( fig_path_and_name , bbox_inches = 'tight' )
print(fig_path_and_name)

show()
