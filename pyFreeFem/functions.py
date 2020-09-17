from pylab import savetxt, array
from numpy import ndarray
from scipy.sparse.linalg import spsolve

import sys
sys.path.append('./../')

from .edpScript import edpScript, edpInput, edpOutput
from .FreeFemTools.edpTools import edp_function, FreeFemize, headerFrame, flagize
from .FreeFemTools.FreeFemStatics import default_variable_names
from .edpScript import *

def InputScript( fespace = None, declare = True, **inputs ) :
    '''
    script = InputScript( fespace = None, declare = True, **inputs )
    '''
    script = edpScript('')

    if fespace is None :
        variable_names = {}
    else :
        variable_names = { '_VhU_' : fespace }

    for input_name, source in inputs.items() :

        if isinstance( source, str ) :
            script += edpInput( name = input_name, data_type = source, variable_names = variable_names, declare = declare )
        else :
            script += edpInput( name = input_name, source = source, variable_names = variable_names, declare = declare )

    return script


def OutputScript( **outputs ) :
    '''
    script = OutputScript( **outputs )
    '''
    script = edpScript('')

    for output_name, output_type in outputs.items() :
        script += edpOutput( name = output_name, data_type = output_type )

    return script

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

    script = InputScript( Th = 'mesh' )

    try :
        u_is_one = u in [ '1', 1 , 1. ]
    except :
        u_is_one = False

    if u is None or u_is_one :
        script += 'Th = ' + edp_function( 'adaptmesh', 'Th', **kwargs )  + ';'
        script_args = dict( Th = Th )

    else :
        script += 'fespace Vh( Th, P1 );'
        script += InputScript( u = 'vector' )
        script += 'Th = ' + edp_function( 'adaptmesh', 'Th', 'u', **kwargs )  + ';'

        script_args = dict( Th = Th, u = u )


    script += OutputScript( Th = 'mesh' )

    _, int_to_label = Th.get_boundary_label_conversion() # Freefem++ can only handle integer boundary labels

    Th = script.get_output( **script_args )['Th']
    Th.rename_boundary(int_to_label)

    return Th

def interpolate( Th, u, P_in, P_out ):
    '''
    w = interpolate( Th, u, P_in, P_out )
    '''

    script = InputScript( Th = Th )

    for fespace, fetype in { 'VhIn' : P_in, 'VhOut' : P_out }.items() :
        script += 'fespace ' + edp_function( fespace, 'Th', fetype ) + ';'

    script += 'VhIn uIn;'
    script += InputScript( uIn = u, declare = False )
    script += 'VhOut uOut = uIn;'
    script += OutputScript( uOut = 'vector' )

    return script.get_output()['uOut']


def get_projector( Th, P_in, P_out ):
    '''
    P = get_projector( Th, P_in, P_out )
    '''

    script = InputScript( Th = Th )

    for fespace, fetype in { 'VhIn' : P_in, 'VhOut' : P_out }.items() :
        script += 'fespace ' + edp_function( fespace, 'Th', fetype ) + ';'

    script += VarfScript(
        GramianOO = 'int2d(Th)( v*u )',
        fespaces = ( 'VhOut', 'VhOut' )
        )

    script += VarfScript(
        GramianIO = 'int2d(Th)( v*u )',
        fespaces = ( 'VhIn', 'VhOut' )
        )

    output = script.get_output()

    GramianOO = output['GramianOO']
    GramianIO = output['GramianIO']

    return spsolve( GramianOO, GramianIO )
