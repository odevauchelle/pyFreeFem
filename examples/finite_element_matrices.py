import sys
sys.path.append('./../')

import pyFreeFem as pyff

script = pyff.edpScript('''
border Circle( t = 0, 2*pi ){ x = cos(t); y = sin(t); }
mesh Th = buildmesh( Circle(20) );
fespace Vh( Th, P1 );
Vh u,v;

fespace Vh2( Th, P2 );
Vh u2, v2;
''' )

# Create and export stiffness matrix
script += pyff.VarfScript( stiffness = 'int2d(Th)( dx(u)*dx(v) +  dy(u)*dy(v) )')

script += pyff.VarfScript( functions = ('u2','v2'), fespaces = ('Vh2', 'Vh2'), stiffness2 = 'int2d(Th)( dx(u2)*dx(v2) +  dy(u2)*dy(v2) )')

stiffness = script.get_output()['stiffness2']


print(stiffness)
