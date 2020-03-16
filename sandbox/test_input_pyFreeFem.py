
import sys
sys.path.append('./../')

import pyFreeFem as pyff

script = pyff.edpScript()
script += pyff.edpInput( name = 'a', data_type = 'int' )
script += pyff.edpInput( name = 'b', data_type = 'real' )
script += 'cout << a*b <<endl;'

print( script.run( a = 2.5, b = 1 ) )
