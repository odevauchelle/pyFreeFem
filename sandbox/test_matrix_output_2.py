
import sys
sys.path.append('./../')

import pyFreeFem as pyff

print(pyff.get_FreeFem_version())
# print(pyff.parse_FreeFem_version('FreeFem++ - version 4.6 (ven. 15 mai 2020 19:54:24 CEST - git no git) 64bits'))

script = pyff.edpScript('''
border Circle( t = 0, 2*pi ){ x = cos(t); y = sin(t); }
mesh Th = buildmesh( Circle(5) );
fespace Vh( Th, P1 );
''')

script += pyff.VarfScript( M = 'int2d(Th)( dx(u)*dx(v) +  dy(u)*dy(v) )' )

print(script.get_edp())

print(script.run())

print( script.get_output()['M'] )
#
# ####################################
#
# with open('output_matrix_FreeFem_4.6.txt', 'r') as the_file :
#     output_freefem = ''.join( the_file.readlines() )
#
# # print(output_freefem)
# #
# print( pyff.FreeFem_str_to_matrix( output_freefem ) )
