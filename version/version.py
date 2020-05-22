import numpy
import matplotlib as mpl
import subprocess

import sys
sys.path.append('./../')

import pyFreeFem as pyff

versions = {
    'pyFreeFem' : pyff.__version__,
    'Python' : sys.version,
    'FreeFem++': pyff.get_FreeFem_version(),
    'Numpy' : numpy.version.version,
    'Matplotlib' : mpl.__version__
    }

version_str = ''

for name, value in versions.items() :

    version_str += '> ' + name + ':\n'
    version_str += str( value )
    version_str += '\n\n'

version_str = version_str[:-2]

with open('version.md','w') as the_file :

    the_file.write('# Versions\n\n')

    the_file.write('Latest test ran with:\n\n')

    the_file.write( "```console\n")
    the_file.write( version_str + '\n' )
    the_file.write( "```\n" )

print( version_str )


################################
#
#   RUN TESTS
#
################################

test_path = './../examples/'

test_files = [
    'quick_example.py',
    'finite_element_matrices.py',
    'mesh_IO.py',
    'boundary_values.py'
    ]

if True :
    for test_file in test_files :

        print('-----------------')
        print('Try to run ' + test_path + test_file  + '?' )

        if input() in ['n', 'no'] :
            break

        subprocess.call( [ 'python3', test_path + test_file ] )
        print('-----------------')

print('All tests ran.')
