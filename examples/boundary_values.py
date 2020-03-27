
import sys
sys.path.append('./../')

from pylab import *
import pyFreeFem as pyff

script = pyff.edpScript('mesh Th = square( 3, 3 );')
script += pyff.OutputScript( Th = 'mesh' )
Th = pyff.adaptmesh( script.get_output()['Th'], hmax = .07 )

u = Th.x*( Th.y - Th.x**3)

fig = plt.figure( constrained_layout = False, figsize = 5*array([1,1.6]) )
gs = fig.add_gridspec( 3, 2 )
ax = fig.add_subplot( gs[:2, :] )
ax_boundary = fig.add_subplot( gs[2:, :], sharex = ax )

ax.tricontourf( Th, u )
Th.plot_triangles( ax = ax, color = 'k', lw = .5, alpha = .2 )


boundaries = Th.get_boundaries()

for label, segments in boundaries.items() :
  print( 'Boundary', label )
  for segment in segments :
    print(*segment)


for boundary_label in 3, 1 :

    boundary = boundaries[boundary_label][0]

    color = ax.plot( Th.x[boundary], Th.y[boundary], clip_on = False, label = boundary_label , lw = 2)[0].get_color()
    ax_boundary.plot( Th.x[boundary], u[boundary], color = color, label = boundary_label )

ax_boundary.set_xlabel('$x$')
ax_boundary.set_ylabel('$u$')
ax.legend( title = 'Boundary', loc = 6 )

ax.axis('equal'); ax.axis('off')
ax.set_yticks([])
#
# fig_path_and_name = './../figures/' + __file__.split('/')[-1].split('.')[0] + '.svg'
# savefig( fig_path_and_name , bbox_inches = 'tight' )
# print(fig_path_and_name)

show()
