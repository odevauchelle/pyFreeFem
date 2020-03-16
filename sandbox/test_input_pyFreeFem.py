
import sys
sys.path.append('./../')

import pyFreeFem as pyff

script = pyff.edpScript()
script += pyff.edpInput( name = 'a', source = 2 )
script += pyff.edpInput( name = 'b', data_type = 'real' )
script += 'a = a*b;'
script += pyff.edpOutput( name = 'a', data_type = 'int' )
script += pyff.edpOutput( name = 'b', data_type = 'real' )

print( script.get_output( b = 1.25 ) )
