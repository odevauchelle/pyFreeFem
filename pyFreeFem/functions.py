from pylab import savetxt, array, shape, diff, real, imag, triu
from numpy import ndarray
from scipy.sparse.linalg import spsolve
from scipy.sparse import lil_matrix


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
        u[0]/1. # u is probably an array
        u_is_one = False

    except :
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

def gradient_matrices( Th, P0 = 'P0', P1 = 'P1' ) :

    '''
    Compute the P1 to P0 gradient matrices, and the areas of the mesh triangles.

    matrices = gradient_matrices( Th, P0 = 'P0', P1 = 'P1' )

    matrices.keys() : ['grad_x', 'grad_y', 'area']

    dxv = matrices['grad_x']*v
    dyv = matrices['grad_y']*v
    '''

    script = InputScript( Th = Th )

    script += edpScript( 'fespace Vh0( Th, ' + P0 + ' ); fespace Vh1( Th, ' + P1 + ' );\n')

    script += VarfScript( grad_x = 'int2d(Th)( u*dx(v) )', grad_y = 'int2d(Th)( u*dy(v) )', fespaces = ('Vh0', 'Vh1') )
    script += VarfScript( area = 'int2d(Th)( u*v )', fespaces = ('Vh0', 'Vh0') )

    matrices = script.get_output()

    matrices['grad_x'] = spsolve( matrices['area'], matrices['grad_x'].T )
    matrices['grad_y'] = spsolve( matrices['area'], matrices['grad_y'].T )

    return matrices


def flux_along_needle( Th, boundary_nodes ) :
    '''
    Compute the ( P0, P0 ) matrix associated to the flux through both sides of a wire.

    q_mat = pyff.flux_along_needle( Th, boundary_nodes )

    Usage:
    q = q_mat*v
    q = q[ boundary_nodes ]
    Q = cumsum(q)
    '''
    matrices = gradient_matrices( Th )

    name_to_index, index_to_name = Th.get_boundary_label_conversion()

    flux = lil_matrix( shape( matrices['grad_x'] ), dtype=np.cfloat ).T

    for i in range( 1, len( boundary_nodes ) ) :

        # identify the two triangles which share edge ( i-1, i )

        edge = [ boundary_nodes[i-1], boundary_nodes[i] ]
        dz = diff( Th.x[edge] + 1j*Th.y[edge] )[0]

        the_sign = -1

        for edge in ( edge, edge[::-1] ) : # two triangles per edge !
            triangle_index = find_triangle_index( Th.triangles, *edge )
            flux[ boundary_nodes[i], triangle_index ] = the_sign*dz
            the_sign *= -1 # different sign for each side of the edge

    return imag( flux ).dot( matrices['grad_x'] ) - real(flux).dot( matrices['grad_y'] )


def integral_along_needle( Th, boundary_nodes ) :
    '''
    Sum P0 contributions along needle.

    cumsum_needle = integral_along_needle( Th, boundary_nodes )
    '''

    cumsum_needle = lil_matrix( triu( [1.]*len(boundary_nodes) ).T )
    boundary_to_Vh = lil_matrix( ( len(Th.x), len(boundary_nodes) ) )

    for j, i in enumerate( boundary_nodes ) :
        boundary_to_Vh[ i, j ] = 1

    return boundary_to_Vh.dot( cumsum_needle.dot( boundary_to_Vh.T ) )


def needle_discharge( Th, boundary_nodes ) :
    '''
    Convenience function to compute the matrix associated to the discharge along a needle.

    Q_mat = needle_discharge( Th, boundary_nodes )
    '''

    return integral_along_needle( Th, boundary_nodes ).dot( flux_along_needle( Th, boundary_nodes ) )
