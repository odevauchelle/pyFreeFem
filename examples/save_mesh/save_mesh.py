
import sys
sys.path.append('./../../')

import pyFreeFem as pyff

script = pyff.edpScript('''
border Circle( t = 0, 2*pi ){ x = cos(t); y = sin(t); }
mesh Th = buildmesh( Circle(5) );
''')

script += pyff.OutputScript( Th = 'mesh' )

Th = script.get_output()['Th']

json_mesh = Th.to_json()
print(json_mesh)

from json import loads

Th = pyff.TriMesh( **loads( json_mesh ) )


Th.plot_triangles( )
Th.plot_boundaries( color = 'red' )

from pylab import *
axis('equal')
show()
