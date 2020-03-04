from pylab import savetxt, array
from numpy import ndarray

import sys
sys.path.append('./../')

from .edpScript import edpScript, edpInput, edpOutput
from .FreeFemTools.edpTools import edp_function

def adaptmesh( Th, u = None, **kwargs ):

    script = edpScript( 'cout << "toto !" << endl;')
    # script += edpInput( name = 'Th', data_type = 'mesh' )

    # if u is None or u is '1' or u is 1 or u is 1. :
    # script += 'Th = ' + edp_function( 'adaptmesh', 'Th', **kwargs )  + ';' # adaptmesh( Th, u, iso = 1 )

    # else :
    #     script += 'fespace Vh( Th, P1 );'
    #     script += edpInput( name = 'u', data_type = 'vector' )
    #     script += 'Th = ' + edp_function( 'adaptmesh', 'Th', 'u', **kwargs )  + ';' # adaptmesh( Th, u, iso = 1 )

    # script += edpOutput( name = 'Th', data_type = 'mesh' )

    print( script.run() )

    # return script.get_output( Th = Th, u = u )['Th']
