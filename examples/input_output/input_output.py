
import sys
sys.path.append('./../')

import pyFreeFem as pyff

script = pyff.InputScript( a = 2, b = 'real' )
script += 'b = a*b;'
script += pyff.InputScript( a = 3, declare = False )
script += pyff.OutputScript( a = 'int', b = 'real' )

print( script.get_output( b = 1.25 ) )
