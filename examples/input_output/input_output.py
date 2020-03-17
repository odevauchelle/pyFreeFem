
import sys
sys.path.append('./../../')

import pyFreeFem as pyff

# script = pyff.InputScript( b = 5.8 )
# print(script.get_edp())
#
#
# script = pyff.InputScript( b = 'real' )
# script += 'cout << b << endl;'
# print( script.run( b = 6.7 ) )
# print( script.run( b = 3 ) )
#
#
# script = pyff.InputScript( a = 2, b = 'real' )
# script += 'b = a*b;'
# script += pyff.InputScript( a = 3, declare = False )
# script += pyff.OutputScript( a = 'int', b = 'real' )
#
# print( script.get_output( b = 1.25 ) )

script = pyff.InputScript( a = 5 )
script += pyff.InputScript( a = 14, declare = False )
script.run()

script = pyff.InputScript( a = 5.1 )
script += pyff.OutputScript( a = 'real' )
print( script.get_output()['a'] )
