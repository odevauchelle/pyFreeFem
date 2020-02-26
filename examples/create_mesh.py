import sys
sys.path.append('./../')

import pyFreeFem as pyff
import matplotlib.pyplot as pp

script = pyff.edpScript('''
border Circle( t = 0, 2*pi ){ x = cos(t); y = sin(t); }
mesh Th = buildmesh( Circle(10) );
''')

script += pyff.edpOutput( type = 'mesh', name = 'Th' )

Th = script.get_output()['Th']

Th.plot_triangles( labels = 'index' )
Th.plot_nodes( labels = 'index', color = 'tab:blue' )
Th.plot_boundaries( color = 'red' )
pp.legend( title = 'Boundary label' )

pp.axis('equal')
pp.axis('off')
pp.xticks([])
pp.yticks([])

pp.savefig( '../figures/' + __file__.split('/')[-1].split('.')[0] + '.svg' , bbox_inches = 'tight' )

pp.show()
