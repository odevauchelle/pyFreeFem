import sys
sys.path.append('./../../')

from pylab import *
import pyFreeFem as pyff

################# PARAMETERS
L = 1.
H = .3
epsilon = .05
npts = 25
n_modes = 30
k = pi/L
levels = logspace( -3, 0, 15 )

def Phi_n( z, n = 1 ) :
    return -cosh( n*k*( 1j*z - H ) )

################ MESH

script = pyff.InputScript( L = L, H = H, epsilon = epsilon, npts = npts )
script +='''
border top( t = L-epsilon, 0 ){ x = t; y = 0; }
border left( t = 0, -(H -epsilon) ){ x = 0; y = t; }
border skirt( t = pi/2, 0 ){ x = epsilon*cos(t); y = -H + epsilon*sin(t); }
border bottom( t = epsilon, L ){x = t; y = -H; }
border right( t = -H, -epsilon ){ x = L; y = t; }
border outlet( t = -pi/2, -pi ){ x = L + epsilon*cos(t); y = epsilon*sin(t); }
mesh Th = buildmesh( top(npts) + left(npts) + skirt(npts) + bottom(npts) + right(npts) + outlet(npts) );
'''
script += pyff.OutputScript( Th = 'mesh' )

Th = script.get_output()['Th']
z_sp = -1j*H # stagnation point
z_out = L

for _ in range(3):
    z = Th.x + 1j*Th.y
    Phi_mesh = exp( n_modes*k*Th.y ) + log( z - z_sp ) + log( z_out - z )
    Th = pyff.adaptmesh( Th, Phi_mesh, iso = 1 )

################ PLOT MESH

# fig = plt.figure(figsize = (6,6))
# gs = fig.add_gridspec( 3, 1 )
# ax_time = fig.add_subplot( gs[2:, :] )
# ax_psi = fig.add_subplot( gs[:2, :], sharex = ax_time)
figure()
ax = gca()
Th.plot_triangles(ax = ax, lw = .7, color = 'k', alpha = .2)
Th.plot_boundaries( color = 'k', clip_on = False )
# ax.tricontourf( Th, real( Phi_mesh ) )
# contours = ax.tricontour( Th, imag( Phi_n( z, n_modes ) ), colors = ['tab:blue'], levels = levels )
xticks([]); yticks([])
ax.axis('equal'); ax.axis('off')

# fig_path_and_name = './../../figures/' + __file__.split('/')[-1].split('.')[0] + '_mesh' + '.svg'
# savefig( fig_path_and_name , bbox_inches = 'tight' )
# print(fig_path_and_name)

################# TRAVEL TIMES
#
# def skirt_time( z, Phi ) :
#     A = real( Phi/z**2 )
#     psi = imag(Phi)
#     phi_in = real(Phi)
#     phi_out = real( -conj( Phi ) ) # symmetry
#     return ( arcsinh( phi_out/psi ) - arcsinh( phi_in/psi ) )/( 4*A )
#
# x_start = []
# t = []
# split = []
#
# z_sp = -1j*H # stagnation point
#
# for contour in contours.collections:
#
#     travel_time_path = 0
#     x_start_path = []
#     split_contour = False
#
#
#     for path in contour.get_paths() :
#         x, y = path.vertices.T
#         z = x + 1j*y
#         phi = real( Phi(z) )
#         ds = sqrt( diff(x)**2 + diff(y)**2 )
#         dt = ds**2/diff(phi)
#         travel_time_path += sum(dt)
#         x_start_path += [ x[0] ]
#
#         if ( abs( z[-1] - z_sp ) - epsilon ) < epsilon :
#             split_contour = True
#             travel_time_path += skirt_time( z[-1] - z_sp, Phi( z[-1] ) - Phi( z_sp ) )
#
#     split += [split_contour]
#     x_start += [ x_start_path[0] ]
#     t += [ travel_time_path ]
#
#
# split = array(split)
# not_split = ~split
# not_split[where(split)[0][-1]] = True
# x_start = array(x_start)
# t = array(t)
#
# figure()
# ax = gca()
# trav_time_color = 'tab:red'
# ax.plot( x_start[not_split], t[not_split], color = trav_time_color, label = 'continuous' )
# ax.plot( x_start[split], t[split], color = trav_time_color, ls = '--', label = 'split' )
# ax.legend(title = 'Contour')
# ax.set_xlabel( 'Starting position' )
# ax.set_ylabel( 'Travel time' )

# fig_path_and_name = './../../figures/' + __file__.split('/')[-1].split('.')[0] + '_trav_time' + '.svg'
# savefig( fig_path_and_name , bbox_inches = 'tight' )
# print(fig_path_and_name)

show()
