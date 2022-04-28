import sys
sys.path.append('./../../')

from pylab import *
import pyFreeFem as pyff

#########################
#
# Create mesh
#
#########################

script = pyff.edpScript('''
border Circle( t = 0, 2*pi ){ x = cos(t); y = sin(t); label = 5; }
mesh Th = buildmesh( Circle(10) );
''')

script += pyff.OutputScript( Th = 'mesh' )

Th = script.get_output()['Th']

#########################
#
# Rename boundary
#
#########################
#
Th.rename_boundary( { 5:'outer wall'} )
print( Th.get_boundary_label_conversion() )


#########################
#
# Return to FreeFem
#
#########################

script = pyff.InputScript( Th = Th )

script += pyff.OutputScript( Th = 'mesh' )

Th = script.get_output()['Th']

#########################
#
# Plot
#
#########################


Th.plot_triangles( lw = 1, color = 'LightGrey')
Th.plot_boundaries( color = 'red' )
legend( title = 'Boundary label' )

axis('equal')
axis('off')
xticks([])
yticks([])

# savefig( './' + __file__.split('/')[-1].split('.')[0] + '_3.svg' , bbox_inches = 'tight' )

show()
