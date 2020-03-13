from pylab import savetxt, array
from numpy import ndarray

import sys
sys.path.append('./../')

from .edpScript import edpScript, edpInput, edpOutput
from .FreeFemTools.edpTools import edp_function, FreeFemize, headerFrame, flagize
from .FreeFemTools.FreeFemStatics import default_variable_names
from .edpScript import *


def VarfBlock( varf, name, functions = None, fespaces = None, FreeFem_name = None, output = True ) :
    '''
    edpBlock = VarfBlock( varf, name, FreeFem_name = None, variable_names = None, output = True )
    '''
    if FreeFem_name is None :
        FreeFem_name = FreeFemize( name, type = 'variable' )

    variable_names = default_variable_names.copy()

    if not functions is None :
        variable_names.update( { '_u_' : functions[0], '_v_' : functions[1]} )

    if not fespaces is None :
        variable_names.update( { '_VhU_' : fespaces[0], '_VhV_' : fespaces[1]} )

    variable_names.update( { '_matrix_name_' : FreeFem_name, '_variational_formulation_' : varf } )

    print(variable_names)


    edp_str = '''
    varf V_matrix_name_( _u_, _v_ ) = _variational_formulation_ ;
    matrix M_matrix_name_ = V_matrix_name_( _VhU_, _VhV_ );
    '''

    for names in variable_names.items() :
        edp_str = edp_str.replace( *names )

    if output :
        output = edpOutput( data_type = 'matrix', name = name, FreeFem_name = FreeFem_name )

    else :
        output = None

    return edpBlock( content = edp_str, output = output )

def VarfScript( functions = None, fespaces = None, **matrices ) :
    '''
    Example :
        edpScript = VarfScript( stiffness = 'int2d(Th)( dx(u)*dx(v) +  dy(u)*dy(v) )', ... )

    The argument's name defines the matrix name in the FreeFem ouput.
    '''

    script = edpScript('')

    for name, varf in matrices.items() :
        script += VarfBlock( name = name, varf = varf, functions = functions, fespaces = fespaces )

    return script

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
        script += 'fespace _Vh_( _Th_, P1 );'
        script += edpInput( name = 'u', data_type = 'vector' )
        script += 'Th = ' + edp_function( 'adaptmesh', 'Th', 'u', **kwargs )  + ';'

        script_args = dict( Th = Th, u = u )


    script += edpOutput( name = 'Th', data_type = 'mesh' )

    _, int_to_label = Th.get_boundary_label_conversion() # Freefem++ can only handle integer boundary labels

    Th = script.get_output( **script_args )['Th']
    Th.rename_boundary(int_to_label)

    return Th
