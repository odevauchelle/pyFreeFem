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
fespace Vh1( Th, P1 );
fespace Vh2( Th, P2 );
'''

script += pyff.VarfScript(
    Gramian1 = 'int2d(Th)( v*u )',
    fespaces = ( 'Vh1', 'Vh1' )
    )

script += pyff.VarfScript(
    Gramian2 = 'int2d(Th)( v*u )',
    fespaces = ( 'Vh2', 'Vh1' )
    )

script += 'Vh2 u2 = x^2 - y^3;'

script += pyff.OutputScript( u2 = 'vector' )

output = script.get_output()

Th = output['Th']
G = None, output['Gramian1'], output['Gramian2']
u2 = output['u2']

proj = spsolve( G[1], G[2] )

print(shape(proj))
print(len( Th.x ))

#################################
#
# Benchmark
#
#################################

script = pyff.InputScript( Th = Th )
script += '''
fespace Vh1( Th, P1 );
fespace Vh2( Th, P2 );
Vh2 u2;
'''
script += pyff.InputScript( u2 = u2, declare = False )

script += 'Vh1 u1 = u2;'

script += pyff.OutputScript( u1 = 'vector' )

u1_ff = script.get_output()['u1']

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

show()
