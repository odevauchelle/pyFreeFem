
import sys
sys.path.append('./../')

import pyFreeFem as pyff
# from pyFreeFem import *

def get_FreeFem_version( parse = True ) :

    version = pyff.run_FreeFem()

    if not parse :
        return version

    else :
        version = version.split('\n')[0]
        version = version.split('version :')[1]
        version = version.strip()

print( get_FreeFem_version() )
