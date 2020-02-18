import sys
sys.path.append('./../')

from pyFreeFem import run_FreeFem
from pyFreeFem.FreeFemTools.edpTools import FreeFemize
from pyFreeFem.FreeFemTools.FreeFemStatics import export_mesh_edp
from pyFreeFem.FreeFemIO import FreeFem_str_to_matrix, parse_FreeFem_output, FreeFem_str_to_mesh

class edpOutput :
    '''
    Defines the type (mesh or matrix) of a FreeFem++ output, its name, and its flag in the FreeFem++ output stream.
    '''
    def __init__( self, type, name, flag = None ) :

        self.type = type
        self.name = name

        if flag is None :
            self.flag = '# ' + FreeFemize( name, type = 'header' )

        else :
            self.flag = flag

    def get_edp( self ) :

        edp = 'cout << "' + self.flag + '" << endl ;'

        if self.type == 'mesh' :
            edp += export_mesh_edp( Th = self.name )

        edp += 'cout << "' + self.flag + '" << endl ;'

        return edp

    def parse( self, FreeFem_output ) :

        if self.type == 'matrix' :
            return FreeFem_str_to_matrix( parse_FreeFem_output( FreeFem_output, self.flag ) )

        elif self.type == 'mesh' :
            return FreeFem_str_to_mesh( parse_FreeFem_output( FreeFem_output, self.flag ) )


class edpBlock :
    '''
    A block of FreeFem++ script. Content in .edp format.
    '''

    def __init__( self, content, name = None, outputs = None, header = None ) :

        self.content = content

        if name is None :
            name = 'unnamed_block'

        self.name = name

        self.outputs = []

        if not outputs is None :
            try :
                self.outputs += outputs
            except :
                self.outputs += [outputs]

        if header is None :
            self.header = FreeFemize( name, type = 'header' )
        else :
            self.header = header

    def get_edp( self ) :

        edp = '/////////////////////////////\n'
        edp += '// ' + self.header + ' START' + '\n'
        edp += '/////////////////////////////\n\n'

        edp += self.content + '\n\n'

        for output in self.outputs :
            edp += output.get_edp()

        edp += '/////////////////////////////\n'
        edp += '// ' + self.header + ' END' + '\n'
        edp += '/////////////////////////////\n\n'


        return edp


class edpScript :
    '''
    Essentially a list of edpBlock objects.
    '''

    def __init__( self, blocks = None ) :

        self.blocks = []

        if not blocks is None :
            self.__add__( blocks )


    def __add__( self, block ) :
        '''
        Argument can be a list of blocks, a string or a block.
        '''
        if type( block ) == type([]) :
            self.blocks += block

        elif type( block ) == type('') :
            self.blocks += [ edpBlock( block ) ]

        else :
            self.blocks += [ block ] # a block


    def get_edp(self) :

        edp = ''
        for block in self.blocks :
            edp += block.get_edp()

        return edp

    def run( self ) :
        return run_FreeFem( self.get_edp() )


    def parse( self, FreeFem_output ) :

        FreeFem_data = {}

        for block in self.blocks :
            for output in block.outputs :
                FreeFem_data[ output.name ] = output.parse( FreeFem_output )

        return FreeFem_data

    def get_FreeFem_data( self ) :
        return self.parse( self.run() )

if __name__ == '__main__' :

    from pylab import *

    mesh_edp = '''
    border Circle( t = 0, 2*pi ){ x = cos(t); y = sin(t); }
    mesh Th = buildmesh( Circle(10) );
    '''

    build_mesh = edpBlock(
        name = 'build_mesh',
        content = mesh_edp,
        outputs = edpOutput( type = 'mesh', name = 'Th' )
        )

    script = edpScript( build_mesh )

    Th = script.get_FreeFem_data()['Th']

    Th.plot_triangles()
    axis('equal')
    show()
