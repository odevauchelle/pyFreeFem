import re
import unicodedata

def FreeFemize( name, type = 'variable' ) :

    name = unicodedata.normalize('NFKD', name )

    name = name.encode('ASCII', 'ignore').decode('utf8')

    name = re.sub('\W+','_', name ) # keeps only alphanumeric or underscore

    if type == 'variable' :
        name = ''.join( map( lambda s: s.capitalize(), name.split('_') ) )

    elif type == 'header' :
        name = ' '.join( ' '.join( name.split('_') ).split() ).upper()

    return name

if __name__ == '__main__' :

    print( FreeFemize( 'toto_ça$$$_vélo$_35' ) )
    print( FreeFemize( 'toto_ça$$$____vélo$_35' , type = 'header' ) )
