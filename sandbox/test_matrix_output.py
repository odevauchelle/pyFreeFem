
import sys
sys.path.append('./../')

import pyFreeFem as pyff

print(pyff.get_FreeFem_version())

script = pyff.edpScript('''
border Circle( t = 0, 2*pi ){ x = cos(t); y = sin(t); }
mesh Th = buildmesh( Circle(5) );
fespace Vh( Th, P1 );
''')

script += pyff.VarfScript( M = 'int2d(Th)( dx(u)*dx(v) +  dy(u)*dy(v) )' )

print(script.get_edp())

print(script.run())

print( script.get_output()['M'] )
