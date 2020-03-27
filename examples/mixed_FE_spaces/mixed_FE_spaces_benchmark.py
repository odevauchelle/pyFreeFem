from pylab import *
from scipy.sparse.linalg import spsolve

import sys
sys.path.append('./../../')

import pyFreeFem as pyff

#################################
#
# Create Mesh
#
#################################

script = pyff.edpScript('mesh Th = square( 8, 8 );')
script += pyff.OutputScript( Th = 'mesh' )

script += '''
fespace Vh( Th, P2 );
Vh u = x^2 - y^3;
'''

script += pyff.OutputScript( u = 'vector' )

output = script.get_output()

Th = output['Th']
u2 = output['u']

proj = pyff.get_projector( Th, 'P2', 'P1' )

print(shape(proj))
print(len( Th.x ))

#################################
#
# Benchmark
#
#################################

u1_ff = pyff.interpolate( Th, u2, 'P2', 'P1' )

figure(figsize = (6,6))

# tricontourf( Th, proj*u2 )
tricontourf( Th, u1_ff )
Th.plot_triangles( color = 'k', lw = .5, alpha = .2 )
tricontour( Th, Th.x**2 - Th.y**3, colors = ['w'], alpha = .3, linestyles = ['--'] )
Th.plot_boundaries( clip_on = False, color = 'k' )

axis('equal')
axis('off')
xticks([]), yticks([])
# legend()

title('Interpolation on P1 space')

# fig_path_and_name = './../../figures/' + __file__.split('/')[-1].split('.')[0] + '.svg'
# savefig( fig_path_and_name , bbox_inches = 'tight' )
# print(fig_path_and_name)

show()
