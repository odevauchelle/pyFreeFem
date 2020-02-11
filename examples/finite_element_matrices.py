import sys
sys.path.append('./../')

import pyFreeFem as pyff

FE_matrix = pyff.stiffness
for key in FE_matrix.keys() :
    print( key + ' : ' + FE_matrix[key] )


edp_str = '''
border Circle( t = 0, 2*pi ){ x = cos(t); y = sin(t); }
mesh Th = buildmesh( Circle(20) );

fespace Vh( Th, P1 );
Vh u,v;
'''

edp_str += pyff.export_matrix_edp( **pyff.stiffness )
FreeFem_output = pyff.run_FreeFem( edp_str )

stiffness_matrix = pyff.FreeFem_str_to_matrix( FreeFem_output, FE_matrix['matrix_name'] )
print(stiffness_matrix)
