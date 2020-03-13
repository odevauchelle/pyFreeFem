from pylab import savetxt, array
from numpy import ndarray

import sys
sys.path.append('./../')

from .edpScript import edpScript, edpInput, edpOutput
from .FreeFemTools.edpTools import edp_function, FreeFemize, headerFrame, flagize
from .FreeFemTools.FreeFemStatics import default_variable_names
from .edpScript import *

def VarfScript( **matrices ) :
    '''
    Example :
        edpScript = VarfScript( stiffness = 'int2d(Th)( dx(u)*dx(v) +  dy(u)*dy(v) )', ... )

    The argument's name defines the matrix name in the FreeFem ouput.
    '''

    script = edpScript('')

    for name, varf in matrices.items() :
        script += VarfBlock( name = name, varf = varf )

    return script


def VarfBlock( varf, name, FreeFem_name = None, variable_names = None, output = True ) :
    '''
    edpBlock = VarfBlock( varf, name, FreeFem_name = None, variable_names = None, output = True )
    '''
    if FreeFem_name is None :
        FreeFem_name = FreeFemize( name, type = 'variable' )

    if variable_names is None :
        variable_names = default_variable_names

    edp_str = '''
    varf Vmatrix_name( base_func, test_func ) = variational_formulation ;
    matrix Mmatrix_name = Vmatrix_name( Vh, Vh ) ;
    '''
    variable_names.update( { 'Mmatrix_name' : FreeFem_name, 'Vmatrix_name' : 'V' + FreeFem_name, 'variational_formulation' : varf } )

    for names in variable_names.items() :
        edp_str = edp_str.replace( *names )

    if output :
        output = edpOutput( data_type = 'matrix', name = name, FreeFem_name = FreeFem_name )

    else :
        output = None

    return edpBlock( content = edp_str, output = output )


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
