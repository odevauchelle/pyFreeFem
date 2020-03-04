from pylab import savetxt, array
from numpy import ndarray

import sys
sys.path.append('./../')

from .edpScript import edpScript, edpInput, edpOutput
from .FreeFemTools.edpTools import edp_function

def adaptmesh( Th, u = None, **kwargs ):
    '''
    TriMesh refinement using FreeFem++'s adaptmesh function.

    Th = adaptmesh( Th, u = None, **kwargs )

    Parameters:
        Th : TriMesh
        u : P1 vector on TriMesh, optional
        other arguments are those of FreeFem++'s adaptmesh function

    Example :

        Th = adaptmesh( Th, hmax = 0.1 )

    '''

    script = edpScript( 'cout << "toto !" << endl;')
    script += edpInput( name = 'Th', data_type = 'mesh' )

    if u is None or u is '1' or u is 1 or u is 1. :
        script += 'Th = ' + edp_function( 'adaptmesh', 'Th', **kwargs )  + ';'

        script_args = dict( Th = Th )


    else :
        script += 'fespace Vh( Th, P1 );'
        script += edpInput( name = 'u', data_type = 'vector' )
        script += 'Th = ' + edp_function( 'adaptmesh', 'Th', 'u', **kwargs )  + ';'

        script_args = dict( Th = Th, u = u )


    script += edpOutput( name = 'Th', data_type = 'mesh' )

    _, int_to_label = Th.get_boundary_label_conversion() # Freefem++ can only handle integer boundary labels

    Th = script.get_output( **script_args )['Th']
    Th.rename_boundary(int_to_label)

    return Th
