import sys
sys.path.append('./../../')

from pylab import *
import pyFreeFem as pyff

################# PARAMETERS
L = 1.
H = .3
epsilon = .05
npts = 15
k = pi/L
levels = logspace( -3, 0, 15 )

def Phi(z) :
    return -cosh( k*( 1j*z - H ) )

################ MESH

script = pyff.InputScript( L = L, H = H, epsilon = epsilon, npts = npts )
script +='''
border top( t = L/2., 0 ){ x = t; y = 0; }
border left( t = 0, -(H -epsilon) ){ x = 0; y = t; }
border skirt( t = pi/2, 0 ){ x = epsilon*cos(t); y = -H + epsilon*sin(t); }
border bottom( t = epsilon, L/2. ){x = t; y = -H; }
border right( t = -H, 0 ){ x = L/2; y = t; }
mesh Th = buildmesh( top(npts) + left(npts) + skirt(npts) + bottom(npts) + right(npts) );
'''
script += pyff.OutputScript( Th = 'mesh' )

Th = script.get_output()['Th']
z = Th.x + 1j*Th.y

################ PLOT MESH

# fig = plt.figure(figsize = (6,6))
# gs = fig.add_gridspec( 3, 1 )
# ax_time = fig.add_subplot( gs[2:, :] )
# ax_psi = fig.add_subplot( gs[:2, :], sharex = ax_time)
figure()
ax = gca()
Th.plot_triangles(ax = ax, lw = .7, color = 'k', alpha = .2)
Th.plot_boundaries( color = 'k', clip_on = False )
contours = ax.tricontour( Th, imag( Phi(z) ), colors = ['tab:blue'], levels = levels )
xticks([]); yticks([])
ax.axis('equal'); ax.axis('off')

# fig_path_and_name = './../../figures/' + __file__.split('/')[-1].split('.')[0] + '_mesh' + '.svg'
# savefig( fig_path_and_name , bbox_inches = 'tight' )
# print(fig_path_and_name)

################# TRAVEL TIMES

def skirt_time( z, Phi ) :
    A = real( Phi/z**2 )
    psi = imag(Phi)
    phi_in = real(Phi)
    phi_out = real( -conj( Phi ) ) # symmetry
    return ( arcsinh( phi_out/psi ) - arcsinh( phi_in/psi ) )/( 4*A )

x_start = []
t = []

z_sp = -1j*H # stagnation point

for contour in contours.collections:

    travel_time_path = 0
    x_start_path = []

    for path in contour.get_paths() :
        x, y = path.vertices.T
        z = x + 1j*y
        phi = real( Phi(z) )
        ds = sqrt( diff(x)**2 + diff(y)**2 )
        dt = ds**2/diff(phi)
        travel_time_path += sum(dt)
        x_start_path += [ x[0] ]

        if ( abs( z[-1] - z_sp ) - epsilon ) < epsilon :
            travel_time_path += skirt_time( z[-1] - z_sp, Phi( z[-1] ) - Phi( z_sp ) )

    x_start += [ x_start_path[0] ]
    t += [ travel_time_path ]

figure()
ax = gca()
trav_time_color = 'tab:red'
ax.plot( x_start, t, color = trav_time_color )
ax.set_xlabel( 'Starting position' )
ax.set_ylabel( 'Travel time' )


show()
