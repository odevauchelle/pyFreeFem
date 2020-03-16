import re
import unicodedata

def input_to_stdin( inputs ) :

    stdin = ''

    for input in inputs :
        stdin += str(input) + ' '

    return stdin

def edp_function( function_name, *args, **kwargs ) :

    edp_str = function_name + '('

    for arg in args :
        edp_str += ' ' + str(arg) + ','

    for arg_name in kwargs.keys() :
        edp_str += ' ' + arg_name + ' = ' + str( kwargs[arg_name] ) + ','

    if edp_str[-1] is ',' :
        edp_str = edp_str[:-1]

    edp_str += ' )'

    return edp_str


def capitalize_first_letter(s):
    try :
        return s[0].upper() + s[1:]
    except :
        return s

def FreeFemize( name, type = 'variable' ) :

    name = unicodedata.normalize('NFKD', name )

    name = name.encode('ASCII', 'ignore').decode('utf8')

    name = re.sub('\W+','_', name ) # keeps only alphanumeric or underscore

    if type == 'variable' :

        if len( name.split('_') ) > 1 :
            name = ''.join( map( lambda s: capitalize_first_letter(s), name.split('_') ) )

    elif type == 'header' :
        name = ' '.join( ' '.join( name.split('_') ).split() ).upper()

    return name

def headerFrame( header ) :
    edp = '\n/////////////////////////////\n'
    edp += '//\n'
    edp += '//    ' + header + '\n'
    edp += '//\n'
    edp += '/////////////////////////////\n\n'
    return edp

def flagize( name ) :
    return '# FLAG > ' + FreeFemize( name, type = 'header' )

if __name__ == '__main__' :

    print( edp_function( 'adaptmesh', 'Th', err = 0.2  )  )

    print( FreeFemize( 'u' ) )
    print( FreeFemize( 'toto_ça$$$_vélo$_35' ) )
    print( headerFrame(FreeFemize( 'toto_ça$$$____vélo$_35' , type = 'header' )) )
    print( flagize( 'début' ) )

    print( input_to_stdin([4.3, 6, 7]) )
