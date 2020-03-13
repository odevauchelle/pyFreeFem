import sys
sys.path.append('./../')

import pyFreeFem as pyff

script = pyff.edpScript('''
border Circle( t = 0, 2*pi ){ x = cos(t); y = sin(t); }
mesh Th = buildmesh( Circle(20) );
fespace Vh( Th, P1 );
Vh u,v;
''' )

# Create and export stiffness matrix
script += pyff.VarfScript( stiffness = 'int2d(Th)( dx(u)*dx(v) +  dy(u)*dy(v) )')

print(script.get_edp())
stiffness = script.get_output()['stiffness']
print(stiffness)
