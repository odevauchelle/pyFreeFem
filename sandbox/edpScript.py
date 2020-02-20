import sys
sys.path.append('./../')

from pyFreeFem import run_FreeFem, TriMesh, savemesh
from pyFreeFem.FreeFemTools.edpTools import FreeFemize, headerFrame, flagize
from pyFreeFem.FreeFemTools.FreeFemStatics import *
from pyFreeFem.FreeFemIO import *

class edpOutput :
    '''
    An output from FreeFem++: mesh, matrix or vector.
    '''
    def __init__( self, type, name, FreeFem_name = None, flag = None, other_variables = None ) :
        '''
        edpOutput( self, type, name, FreeFem_name = None, flag = None, other_variables = None )

        Arguments:
            type (str): 'mesh', 'matrix' or 'vector'
            name (str): name under which the output will be returned
            FreeFem_name (str): name of the variable in the edp file
            flag (str): flag signalling the output in the FreeFem++ stream
            other_variables (dict) : translation dictionnary for FreeFem++ variables
        '''

        self.type = type
        self.name = name

        if FreeFem_name is None :
            self.FreeFem_name = FreeFemize( name, type = 'variable' )
        else :
            self.FreeFem_name = FreeFem_name

        if flag is None :
            self.flag = flagize( name )
        else :
            self.flag = flag

        if other_variables is None :
            self.other_variables = {}
        else :
            self.other_variables = other_variables

    def get_edp( self ) :

        if self.type == 'matrix' :
            edp = export_matrix_edp( create_varf = False, create_and_add_flags = False, Mmatrix_name = self.FreeFem_name, **self.other_variables )

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

    def __init__( self, source, name, FreeFem_name = None, source_type = None, declare_variable = True ) :
        '''
        edpInput( source, name, FreeFem_name = None, source_type = None )

        Arguments:
            source ( TriMesh, list, array, sparse matrix ) : to be loaded into FreeFem++
            name (str) : input name
            FreeFem_name (str) : name of the corresponding FreeFem++ variable
            source_type (str) : mesh, vector, matrix or number
            declare_variable (bool) : whether to declare the FreeFem++ variable or not
        '''

        self.source = source
        self.name = name
        self.declare_variable = declare_variable

        if source_type is None :

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
            self.type = source_type


        if FreeFem_name is None :
            self.FreeFem_name = FreeFemize( name, type = 'variable' )
        else :
            self.FreeFem_name = FreeFem_name

    def get_edp( self, **kwargs ) :

        edp_str = ''

        if self.type is 'mesh' :

            temp_mesh_file = NamedTemporaryFile( suffix = '.msh' )

            savemesh( filename = temp_mesh_file.name, mesh = self.source )

            substitutions = { 'mesh_file_name' : temp_mesh_file.name, 'Th' : self.FreeFem_name }

            if self.declare_variable :
                edp_str += 'mesh Th;\n'

            edp_str += '''
            Th = readmesh( "mesh_file_name" ) ;
            '''

            for key in substitutions.keys() :
                edp_str = edp_str.replace( key, substitutions[key] )

        return edp_str



class edpBlock :
    '''
    A block of FreeFem++ script. Content in .edp format.
    '''

    def __init__( self, content = None, name = None, output = None, header = None ) :

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

        self.outputs = []

        if not output is None :
            try :
                self.outputs += output
            except :
                self.outputs += [output]

        if header is None :
            self.header = FreeFemize( name, type = 'header' )
        else :
            self.header = header

    def get_edp( self ) :

        edp = headerFrame( self.header + ' START' )

        edp += self.content + '\n\n'

        for output in self.outputs :
            edp += output.get_edp()

        edp += headerFrame( self.header + ' END' )

        return edp


class edpScript :
    '''
    Essentially a list of edpBlock objects.
    '''

    def __init__( self, name = None, blocks = None ) :

        '''
        edpScript( blocks = None ) :

        Arguments:
            blocks (list, string or edpBlock)

        '''

        self.name = name
        self.blocks = []

        if not blocks is None :
            self.__add__( blocks )


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
            self.blocks += block

        elif type( other ) == type('') :
            self.blocks += [ edpBlock( block ) ]

        elif type( other ) == edpBlock :
            self.blocks += [ other ] # a block

        elif type( other ) == edpOutput :
            self.blocks += [ edpBlock( output = other ) ] # a block

    def __add__( self, other ) :

        self.add( other )

        return self

    def get_edp(self) :

        edp = ''

        if not self.name is None :
            edp += headerFrame( FreeFemize( self.name, type = 'header' ) )

        for block in self.blocks :
            edp += block.get_edp()

        return edp

    def run( self, **kwargs ) :
        return run_FreeFem( self.get_edp(), **kwargs )


    def parse( self, FreeFem_output ) :

        FreeFem_data = {}

        for block in self.blocks :
            for output in block.outputs :
                FreeFem_data[ output.name ] = output.parse( FreeFem_output )

        return FreeFem_data

    def get_output( self ) :
        return self.parse( self.run() )

if __name__ == '__main__' :

    from pylab import *

    script = edpScript()

    mesh_edp = '''
    border Circle( t = 0, 2*pi ){ x = cos(t); y = sin(t); }
    mesh Th = buildmesh( Circle(10) );
    fespace Vh( Th, P1 );
    Vh u,v;
    '''

    script += edpBlock(
        name = 'build_mesh',
        content = mesh_edp,
        )

    script += edpOutput( type = 'mesh', name = 'Th' )

    script += edpBlock(
        name = 'adapt_mesh',
        content = 'Th = adaptmesh( Th, 1, hmax = .2 );',
        output = edpOutput( type = 'mesh', name = 'Th2', FreeFem_name = 'Th' )
        )

    script += edpBlock(
        name = 'vector output',
        content = 'u = x^2;',
        output = edpOutput( type = 'vector', name = 'x2', FreeFem_name = 'u' )
        )

    script += edpBlock(
        name = 'stiffness_matrix',
        content = create_varf_matrix( **stiffness ),
        output = edpOutput( type = 'matrix', name = 'stiffness', FreeFem_name = 'Mstiffness', other_variables = stiffness )
        )

    # script = edpScript( name = 'This is a title' ) + script

    # print(script.get_edp())
    # print(script.run())

    FFdata = script.get_output()
    #

    input = edpInput( name = 'Th', source = FFdata['Th2'] )
    print( input.get_edp() )


    # FFdata['Th'].plot_triangles()
    tricontourf( FFdata['Th2'], FFdata['x2'] )

    FFdata['Th2'].plot_triangles( alpha = .2 )

    axis('equal')
    show()
