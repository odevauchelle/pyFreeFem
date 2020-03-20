import sys
sys.path.append('./../../')

from pylab import *
import pyFreeFem as pyff
from matplotlib.tri import LinearTriInterpolator
from numpy.ma import is_masked

################# PARAMETERS
L = 1.
H = .3
epsilon = .05
npts = 20
n_modes = 500
k = pi/L
levels = logspace( -3, 0, 20 )

def Phi_n( z, n = 1 ) :

    shallow_z = imag(z) > -8./( n*k ) # deeper than this, the exponential virtually vanishes

    try : # z is an array
        Phi = z*0
        Phi[shallow_z] = -cosh( n*k*( 1j*z[shallow_z] - H ) )
        return Phi

    except : # z is a single complex number
        if shallow_z :
            return -cosh( n*k*( 1j*z - H ) )
        else :
            return 0


def Phi_sum( z ) :

    Phi = 0.*z

    for n in range( 1, n_modes ) :
        Phi += 2*Phi_n( z, n )/( n*k*sinh( n*k*H ) )

    return Phi

################ MESH

script = pyff.InputScript( L = L, H = H, epsilon = epsilon, npts = npts )
script +='''
border top( t = L, epsilon ){ x = t; y = 0; }
border outlet( t = 0, -pi/2 ){ x = epsilon*cos(t); y = epsilon*sin(t); }
border left( t = -epsilon, -( H -epsilon) ){ x = 0; y = t; }
border skirt( t = pi/2, 0 ){ x = epsilon*cos(t); y = -H + epsilon*sin(t); }
border bottom( t = epsilon, L - epsilon ){x = t; y = -H; }
border skirtRight( t = pi, pi/2 ){ x = L + epsilon*cos(t); y = -H + epsilon*sin(t); }
border right( t = -H+epsilon, 0 ){ x = L; y = t; }
mesh Th = buildmesh( top(npts) + outlet(npts) + left(npts) + skirt(npts) + bottom(npts) + skirtRight(npts) + right(npts) );
'''
script += pyff.OutputScript( Th = 'mesh' )
Th = script.get_output()['Th']

for _ in range(3):
    Th = pyff.adaptmesh( Th, iso = 1, hmax = epsilon/2, err = 5e-3 )

################ Phi Interpolator

Phi_value = Phi_sum( Th.x + 1j*Th.y )

phi_int, psi_int = LinearTriInterpolator( Th, real( Phi_value ) ), LinearTriInterpolator( Th, imag( Phi_value ) )

def Phi( z ) :

    x, y = real(z), imag(z)
    Phi_int = phi_int( x, y ) + 1j*psi_int( x, y )

    if is_masked(Phi_int) : # the interpolator fails
        return Phi_sum( z )

    else :
        return Phi_int

################ PLOT MESH

# fig = plt.figure(figsize = (6,6))
# gs = fig.add_gridspec( 3, 1 )
# ax_time = fig.add_subplot( gs[2:, :] )
# ax_psi = fig.add_subplot( gs[:2, :], sharex = ax_time)
figure(figsize = (7,3.))
ax = gca()
Th.plot_triangles(ax = ax, lw = .7, color = 'k', alpha = .2)
Th.plot_boundaries( color = 'k', clip_on = False )
contours = ax.tricontour( Th, imag( Phi_value ), colors = ['tab:blue'], levels = levels )
xticks([]); yticks([])
ax.axis('equal'); ax.axis('off')
#
ax.set_title( str( n_modes ) + ' modes' )
fig_path_and_name = './../../figures/' + __file__.split('/')[-1].split('.')[0] + '_mesh' + '.svg'
savefig( fig_path_and_name , bbox_inches = 'tight' )
print(fig_path_and_name)

############### TRAVEL TIMES

z_sp_left = -1j*H # stagnation point
z_sp_right = L - 1j*H # stagnation point
z_out = 0

def skirt_time( z, Phi ) :
    A = real( Phi/z**2 )
    psi = imag( Phi )
    phi_in = real( Phi )
    phi_out = real( -conj( Phi ) ) # symmetry
    return abs( ( arcsinh( phi_out/psi ) - arcsinh( phi_in/psi ) )/( 4*A ) )

x_start = []
t = []
is_split = []

for contour in contours.collections[:-1]:

    travel_time_path = 0.
    x_start_path = []
    is_split_contour = False

    for path in contour.get_paths() :

        x, y = path.vertices.T
        z = x + 1j*y
        phi = real( Phi( z ) )
        ds = sqrt( diff(x)**2 + diff(y)**2 )
        dt = ds**2/diff(phi)
        travel_time_path += sum(dt)

        for z_sp in ( z_sp_left, z_sp_right ) :

            z_in = z[0]

            if ( abs( z_in - z_sp ) - epsilon ) < epsilon**2 :
                travel_time_path += skirt_time( z_in - z_sp, Phi( z_in ) -  Phi( z_sp ) )
                is_split_contour = True

    x_start += [ contour.get_paths()[-1].vertices[-1][0] ]
    t += [ travel_time_path ]
    is_split += [ is_split_contour ]


is_split = array(is_split)
not_is_split = ~is_split
not_is_split[where(is_split)[0][-1]] = True
x_start = array(x_start)
t = array(t)

figure()
ax = gca()
trav_time_color = 'tab:red'
ax.plot( x_start[not_is_split], t[not_is_split], color = trav_time_color, label = 'continuous' )
ax.plot( x_start[is_split], t[is_split], color = trav_time_color, ls = '--', label = 'is_split' )
ax.legend(title = 'Contour')
ax.set_xlabel( 'Starting position' )
ax.set_ylabel( 'Travel time' )

# fig_path_and_name = './../../figures/' + __file__.split('/')[-1].split('.')[0] + '_trav_time' + '.svg'
# savefig( fig_path_and_name , bbox_inches = 'tight' )
# print(fig_path_and_name)

show()
