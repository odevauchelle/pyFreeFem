
import sys
sys.path.append('./../')

import pyFreeFem as pyff
from pylab import array
from numpy import float64

A = array([[2.0]]*2)

a = A[0,0]
print(type(a) is float64)
print(type(2.6) is float)

script = pyff.InputScript( a = a )
script += pyff.OutputScript( a = 'real')

print( script.get_output() )
