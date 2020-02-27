from pylab import savetxt, array
from numpy import ndarray

import sys
sys.path.append('./../')

from .TriMesh import TriMesh
from .FreeFemTools.edpTools import FreeFemize, headerFrame, flagize
from .FreeFemTools.FreeFemStatics import *
from .FreeFemIO import *

default_variable_names = {
    'base_func' : 'u',
    'test_func' : 'v'
    }

class edpOutput :
    '''
    An output from FreeFem++: mesh, matrix or vector.
    '''
    def __init__( self, data_type, name, FreeFem_name = None, flag = None, variable_names = None ) :
        '''
        edpOutput( self, type, name, FreeFem_name = None, flag = None, variable_names = None )

        Arguments:
            type (str): 'mesh', 'matrix' or 'vector'
            name (str): name under which the output will be returned
            FreeFem_name (str): name of the variable in the edp file
            flag (str): flag signalling the output in the FreeFem++ stream
            variable_names (dict) : translation dictionnary for FreeFem++ variables
        '''

        self.type = data_type
        self.name = name

        if FreeFem_name is None :
            self.FreeFem_name = FreeFemize( name, type = 'variable' )
        else :
            self.FreeFem_name = FreeFem_name

        if flag is None :
            self.flag = flagize( name )
        else :
            self.flag = flag

        if variable_names is None :
            self.variable_names = default_variable_names
        else :
            self.variable_names = variable_names

    def get_edp( self ) :

        if self.type == 'matrix' :
            variable_names = self.variable_names
            variable_names.update( { 'Mmatrix_name' : self.FreeFem_name } )
            edp = export_matrix_edp( create_varf = False, create_and_add_flags = False, **variable_names )

        elif self.type == 'vector' :
            edp = export_vector_edp( base_func = self.FreeFem_name )

        elif self.type == 'mesh' :
            edp = export_mesh_edp( Th = self.FreeFem_name )

        return add_flags( edp, self.flag )

    def parse( self, FreeFem_output ) :

        if self.type == 'matrix' :
            return FreeFem_str_to_matrix( parse_FreeFem_output( FreeFem_output, self.flag  ) )

        elif self.type == 'vector' :
            return FreeFem_str_to_vector( parse_FreeFem_output( FreeFem_output, self.flag ) )

        elif self.type == 'mesh' :
            return FreeFem_str_to_mesh( parse_FreeFem_output( FreeFem_output, self.flag ) )

class edpInput :
    '''
    Input into FreeFem.
    '''

    def __init__( self, name, source = None, FreeFem_name = None, data_type = None, tempfile = None, declare = True, variable_names = None ) :
        '''
        edpInput( source, name, FreeFem_name = None, data_type = None )

        Arguments:
            source ( TriMesh, list, array, sparse matrix ) : to be loaded into FreeFem++
            name (str) : input name
            FreeFem_name (str) : name of the corresponding FreeFem++ variable
            data_type (str) : mesh, vector, matrix or number
        '''

        self.source = source
        self.name = name
        self.tempfile = tempfile # including when it's None
        self.declare = declare

        if variable_names is None :
            self.variable_names = default_variable_names
        else :
            self.variable_names = variable_names

        if FreeFem_name is None :
            self.FreeFem_name = FreeFemize( name, type = 'variable' )
        else :
            self.FreeFem_name = FreeFem_name

        if data_type is None :

            if type(source) is TriMesh :
                self.type = 'mesh'

            elif type(source) is int :
                self.type = 'int'

            elif type(source) is float :
                self.type = 'real'

            elif type(source) in [ list, array, ndarray ] :
                self.type = 'vector'

            elif type(source) is default_sparse_matrix :
                self.type = 'matrix'

            else :
                self.type = None

        else :
            self.type = data_type

    def get_edp( self, **kwargs ) :

        edp_str = ''

        if self.source is None :
            source = kwargs[ self.name ] # the source is input when calling get_edp

        else :
            source = self.source

        if self.type is 'mesh' :

            if self.tempfile is None :
                self.tempfile = NamedTemporaryFile( suffix = '.msh' )

            savemesh( filename = self.tempfile.name, mesh = source )

            if self.declare :
                edp_str += 'mesh Th;\n'

            edp_str += 'Th = readmesh( "mesh_file_name" ) ;\n'

            variable_names = self.variable_names
            variable_names.update( { 'mesh_file_name' : self.tempfile.name, 'Th' : self.FreeFem_name } )

            for key in variable_names.keys() :
                edp_str = edp_str.replace( key, variable_names[key] )

        elif self.type is 'vector' :

            if self.tempfile is None :
                self.tempfile = NamedTemporaryFile( suffix = '.ffv' )

            savetxt( self.tempfile.name, source ) # using the file handle would be better, but then writing doesn't complete

            if self.declare :
                edp_str += 'Vh vector_name;\n'

            edp_str += '''
                {
                    ifstream InputFile("vector_file_name");

                    for(int i = 0; i < vector_name.n; i++)
                        {
                        InputFile >> vector_name[][i];
                        }
                }
            '''

            variable_names = self.variable_names
            variable_names.update( { 'vector_file_name' : self.tempfile.name, 'vector_name' : self.FreeFem_name } )

            for key in variable_names.keys() :
                edp_str = edp_str.replace( key, variable_names[key] )

        return edp_str

class edpBlock :

    '''
    A block of FreeFem++ script. Content in .edp format.
    '''

    def __init__( self, content = None, name = None, input = None, output = None, header = None ) :

        '''
        edpBlock(content, name = None, outputs = None, header = None ) :

        Arguments :
            content (str) : edp script for FreeFem++
            name (str) : name of the block
            output (list or edpOutput) : output to be expected from this block
            header (str) : header of the block in the FreeFem++ script
        '''

        if content is None :
            self.content = ''
        else :
            self.content = content

        if name is None :
            name = 'unnamed_block'

        self.name = name

        self.input = []
        if not input is None :
            try :
                self.input += input
            except :
                self.input += [input]

        self.output = []
        if not output is None :
            try :
                self.output += output
            except :
                self.output += [output]

        if header is None :
            self.header = FreeFemize( name, type = 'header' )
        else :
            self.header = header

    def get_edp( self, **kwargs_input ) :

        edp = headerFrame( self.header + ' START' )

        for input in self.input :
            edp += input.get_edp( **kwargs_input )

        edp += self.content + '\n\n'

        for output in self.output :
            edp += output.get_edp()

        edp += headerFrame( self.header + ' END' )

        return edp


class edpScript :
    '''
    Essentially a list of edpBlock objects.
    '''

    def __init__( self, blocks = None, name = None ) :

        '''
        edpScript( blocks = None ) :

        Arguments:
            blocks (list, string or edpBlock)

        '''

        self.blocks = []

        if not blocks is None :
            self.add( blocks )

        self.name = name

    def add( self, other ) :
        '''
        Argument can be:
            Another edpScript
            An edpBlock
            A list of edpBlocks
            An edpOutput
        '''

        if type( other ) == edpScript :
            # name property of other is lost
            self.blocks += other.blocks

        if type( other ) == type([]) :
            self.blocks += other

        elif type( other ) == type('') :
            self.blocks += [ edpBlock( other ) ]

        elif type( other ) == edpBlock :
            self.blocks += [ other ] # a block

        elif type( other ) == edpOutput :
            self.blocks += [ edpBlock( output = other ) ] # a block

        elif type( other ) == edpInput :
            self.blocks += [ edpBlock( input = other ) ] # a block

    def __add__( self, other ) :

        self.add( other )

        return self

    def get_edp( self, **kwargs_input ) :

        edp = ''

        if not self.name is None :
            edp += headerFrame( FreeFemize( self.name, type = 'header' ) )

        for block in self.blocks :
            edp += block.get_edp( **kwargs_input )

        return edp

    def clean_temp_files(self, verbose = False) :

        for block in self.blocks :

            for input in block.input :

                if not input.tempfile is None :

                    if verbose :
                        print( 'Erasing temporary file ' + input.tempfile.name )

                    input.tempfile.close()

    def run( self, **kwargs_input ) :
        freefem_output = run_FreeFem( self.get_edp( **kwargs_input ) )
        self.clean_temp_files()
        return freefem_output

    def parse( self, FreeFem_output ) :

        FreeFem_data = {}

        for block in self.blocks :
            for output in block.output :
                FreeFem_data[ output.name ] = output.parse( FreeFem_output )

        return FreeFem_data

    def get_output( self, **kwargs_input ) :
        return self.parse( self.run( **kwargs_input ) )


def VarfBlock( varf, name, FreeFem_name = None, variable_names = None, output = True ) :

        if FreeFem_name is None :
            FreeFem_name = FreeFemize( name, type = 'variable' )

        if variable_names is None :
            variable_names = default_variable_names

        edp_str = '''
        varf Vmatrix_name( base_func, test_func ) = variational_formulation ;
        matrix Mmatrix_name = Vmatrix_name( Vh, Vh ) ;
        '''
        variable_names.update( { 'Mmatrix_name' : FreeFem_name, 'Vmatrix_name' : 'V' + FreeFem_name, 'variational_formulation' : varf } )

        for key in variable_names.keys() :
            edp_str = edp_str.replace( key, variable_names[key] )

        if output :
            output = edpOutput( data_type = 'matrix', name = name, FreeFem_name = FreeFem_name )

        else :
            output = None

        return edpBlock( content = edp_str, output = output )




# if __name__ == '__main__' :
#
#     from pylab import *
#
#     script = edpScript( name = 'build_mesh',
#         blocks = edpBlock(
#             content = '''
#             border Circle( t = 0, 2*pi ){ x = cos(t); y = sin(t); }
#             mesh Th = buildmesh( Circle(15) );
#             ''')
#         )
#
#     script += edpBlock(
#         name = 'adapt_mesh',
#         content = 'Th = adaptmesh( Th, 1, hmax = .01 );',
#         output = edpOutput( type = 'mesh', name = 'Th_refined', FreeFem_name = 'Th' )
#         )
#
#     Th = script.get_output()['Th_refined']
#
#     script = edpScript('')
#
#     u = sin( Th.x*pi*5 + 2*Th.y**3 )
#
#     script += edpBlock( input = edpInput( name = 'Th', source = Th ) )
#
#     script += '''
#     fespace Vh( Th, P1 );
#     '''
#
#     script += edpBlock( input = edpInput( name = 'u', FreeFem_name = 'u', source = u ) )
#
#     script += edpBlock(
#         name = 'vector output',
#         output = edpOutput( type = 'vector', name = 'u', FreeFem_name = 'u' )
#         )
#
#
#
#     # script += edpBlock(
#     #     name = 'stiffness_matrix',
#     #     content = create_varf_matrix( **stiffness ),
#     #     output = edpOutput( type = 'matrix', name = 'stiffness', FreeFem_name = 'Mstiffness', variable_names = stiffness )
#     #     )
#
#     FFdata = script.get_output()
#     tricontourf( Th, FFdata['u'] )
#     Th.plot_triangles( alpha = .2, color = 'w', lw = 1 )
#
#     axis('equal')
#     show()
