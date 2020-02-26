import sys
sys.path.append('./../')

import pyFreeFem as pyff

script = pyff.edpScript( 'cout << "Hello world!" << endl;' )

print( script.run() )
