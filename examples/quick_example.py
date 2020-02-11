import sys
sys.path.append('./../')

import pyFreeFem as pyff

FreeFem_output = pyff.run_FreeFem( 'cout << "Hello world!" << endl;' )

print(FreeFem_output)
