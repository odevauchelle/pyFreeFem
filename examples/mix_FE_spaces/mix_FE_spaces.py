from pylab import *

import sys
sys.path.append('./../../')

import pyFreeFem as pyff

#################################
#
# Create Mesh
#
#################################

script = pyff.edpScript('mesh Th = square(20, 20);')
script += pyff.edpOutput( data_type = 'mesh', name = 'Th' )
Th = script.get_output()['Th']

#################################
#
# Get matrices
#
#################################

script = pyff.edpScript()
script += pyff.edpInput( 'Th', Th )

script += '''
fespace Vh( Th, P2 );
Vh u, v;
'''

script += pyff.VarfScript(
    stiffness = 'int2d(Th)( dx(u)*dx(v) +  dy(u)*dy(v) )',
    Grammian = 'int2d(Th)( u*v )',
    boundary_Grammian = 'int1d( Th, 1, 2, 3, 4 )( u*v )'
    )

M = script.get_output()

Th.plot_triangles(color = 'k', lw = .5, alpha = .2)
Th.plot_boundaries()

axis('equal')
axis('off')
xticks([]), yticks([])
legend()

show()
