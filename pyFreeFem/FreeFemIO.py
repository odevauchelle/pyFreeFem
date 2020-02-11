#!/usr/bin/env python
# -*- coding: utf-8 -*-

# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <https://www.gnu.org/licenses/>.
#
# Olivier Devauchelle
#
# Many thanks to Eric Lajeunesse, Anaïs Abramian, Valentin Jules & Hugo Chauvet.
#
# When used in a scientific publication, please cite:
#
# Boltzmann Distribution of Sediment Transport, A. Abramian, O. Devauchelle,
# G. Seizilles, E. Lajeunesse, Physical Review Letters, 123, 014501, 2019

import numpy as np
import re
import csv
import subprocess
from scipy.sparse import csr_matrix
from scipy.sparse.linalg import spsolve

from .meshing import TriMesh, Boundary
from .FreeFemStatics import *

def FreeFem_str_to_matrix( FreeFem_str, matrix_name, raw = False ) :
    '''
    '''

    flag = '# MATRIX ' + matrix_name + '\n'

    FreeFem_lines = FreeFem_str.split( flag )[1].split('\n')
    nb_row, nb_col, is_symmetric, nb_coef = map( lambda x : int(x), FreeFem_lines[3].split() )
    I, J, coef = loadstr( '\n'.join(  FreeFem_lines[4:-1] ) ).T
    I = np.array( list( map( lambda x: int(x), I ) ) ) - 1
    J = np.array( list( map( lambda x: int(x), J ) ) ) - 1
    coef = np.array( list( map( lambda x: float(x), coef ) ) )

    if raw :
        return ( coef, (I, J) ), ( nb_row, nb_col )
    else :
        return csr_matrix( ( coef, (I, J) ), ( nb_row, nb_col ) )

def FreeFem_str_to_mesh( FreeFem_str, simple_boundaries = True ) :

    '''
    '''

    FreeFem_str = '\n' + FreeFem_str + '\n'

    flag = {
        'nodes' : '\n# NODES\n',
        'boundaries' : '\n# BOUNDARIES\n',
        'triangles' : '\n# TRIANGLES\n'
        }

    FreeFem_mesh = {}

    for key in ['triangles', 'boundaries'] :
        FreeFem_mesh[key] = loadstr( FreeFem_str.split( flag[key] )[1], dtype = 'int' )

    for key in ['nodes'] :
        FreeFem_mesh[key] = loadstr( FreeFem_str.split( flag[key] )[1], dtype = 'float' )

    x, y, node_labels = FreeFem_mesh['nodes'].T
    node_labels = list( map( lambda x: int(x), node_labels ) )

    triangles = FreeFem_mesh['triangles'][:,:-1]
    triangle_labels = FreeFem_mesh['triangles'][:,-1]

    boundaries = FreeFem_to_boundaries( FreeFem_mesh['boundaries'] )

    mesh = TriMesh(
        x, y,
        triangles = triangles,
        node_labels = node_labels,
        triangle_labels = triangle_labels,
        boundaries = boundaries
        )

    return mesh

######################

def find_disjunctions( seg_list ) :
    '''
    Find disjunctions (holes) in a list of segments.

    Parameters :
    ----------
    seg_list : numpy.array of segments
        list of couples of nodes coordinates.

    Returns :
    ----------
    disjuctions : numpy.array of indices
        Indices of the disjunction segments

    Example :
    ----------
    disjuctions = find_disjunction( array( [ [1,12], [12,3], [3,7], [4, 14], [14,1] ] ) )
    >>> [3]

    '''

    if seg_list is [] :
        return []

    else :
        return ( seg_list[ 1:-1, 0 ] - seg_list[ 0:-2, 1 ] ).nonzero()[0] + 1

def stitch_segment_list( seg_list, n ) :

    '''
    Cut a segment list at location 'n', finds the next occurence of the end node, and stich the rest of the list there.

    Parameters :
    ----------
    seg_list : numpy.array of segments
        List of couples of nodes coordinates.
    n : integer
        Location of the disjuction to be mended.

    Returns :
    ----------
    mended_segment_list : numpy.array of indices
        Indices of the disjunction segments
    mending : Boolean
        Whether the list of segment was mended

    Example :
    ----------
    stitch_segment_list( array( [ [1,12], [12,3], [7, 14], [3, 7], [14, 1 ] ] ), 2 )
    >>> (array([[ 1, 12], [12,  3], [ 3,  7], [14,  1], [ 7, 14]]), True)
    '''

    try :
        n2 = n + np.where( seg_list [ n:, 0  ] == seg_list[ ( n - 1 ), 1 ] )[0][0]
        return np.concatenate( ( seg_list[ 0 : n ], seg_list[ n2 : ], seg_list[ n : n2 ] )  ), True

    except :
        print('No possible reconnection of segment list. Returning initial segment list.')
        return seg_list, False

def reorder_boundary( seg_list ) :
    '''
    Reorder a list of segments until there is no disjunction left, or until the remaining disjunctions cannot be solved.

    Parameters :
    ----------
    seg_list : numpy.array of segments
        List of couples of nodes coordinates.

    Returns :
    ----------
    ordered_segment_list : numpy.array of indices
        List of ordered segments, with the least possible disjuctions

    Example :
    ----------
    reorder_boundary( array( [ [1,12], [12,3], [7, 14], [3, 7], [14, 1 ] ] ) )
    >>> array( [[ 1 12], [12  3], [ 3  7], [ 7 14], [14  1]] )
    '''

    keep_trying = True

    while keep_trying :

        keep_trying = False

        for disjunction_index in find_disjunctions( seg_list ) :

            seg_list, mending = stitch_segment_list( seg_list, disjunction_index )

            if mending :
                keep_trying = True
                break

    return seg_list

def FreeFem_to_boundaries( seg_list, reorder = True ) :
    '''
    Turns a list of labelled segments into a list of boundaries.
    '''

    boundaries = []

    for label in set( seg_list[:,2] ) :

        segments = seg_list[ np.where( seg_list[:,2] == label ), :2 ]

        if reorder :
            segments = reorder_boundary( segments )

        boundaries += [ Boundary( segments = segments, label = label ) ]

    return boundaries

def loadstr( data_str, delimiter = None, dtype = 'float', skip_rows = 0 ) :

    '''
    An approximation of np.loadtxt that applies to str.
    '''


    if dtype is 'float' :
        def conversion_try( x ) :
            return float( x )

    elif dtype is 'int' :
        def conversion_try( x ) :
            return int( x )

    else :
        def conversion_try( x ) :
            return x

    def conversion(x) :
        try :
            return conversion_try(x)
        except :
            return x

    data = []

    for data_line in data_str.split('\n')[skip_rows:] :
        data += [ list( map( conversion, data_line.split( delimiter ) ) ) ]

    return np.array(data)

def run_FreeFem( edp_str ) :
    '''
    Run FreeFem++ on script edp_str, and returns Popen output.
    '''

    edp_str = edp_str.replace( '"', "\"'\"" )

    command = 'FreeFem++ -v 0 <( printf "' + edp_str + '" )'

    # output = subprocess.check_output( [command], shell = True, stderr=subprocess.STDOUT, executable="/bin/bash", bufsize = 4096 )
    #
    # return output.decode('utf-8')
    #
    # subprocess.run( ['stdbuf', '-o0', '&&'], shell = True )
    # output = subprocess.run( [command], stdout=subprocess.PIPE, shell = True, executable="/bin/bash" )
    # return output.stdout.decode('utf-8')

    #
    # proc = Popen(['stdbuf', '-o0']

    with subprocess.Popen( [command], stdout = subprocess.PIPE, shell = True, executable="/bin/bash", bufsize = 0 ) as proc :

        output, error = proc.communicate()

        return output.decode('utf-8')
        proc.wait()
        proc.stdout.flush()
        # proc.kill()



if __name__ == '__main__' :

    from pylab import *

    edp_str = '''
    real smallRadius = .3;
    border outerCircle( t = 0, 2*pi ){ x = cos(t); y = 0.8*sin(t); }
    border innerCircle( t = 2*pi, 0 ){ x = .5 + smallRadius*cos(t); y = smallRadius*sin(t); }
    mesh Th = buildmesh( outerCircle(100) + innerCircle(40) );

    fespace Vh( Th, P1 );
    Vh u,v;
    '''

    edp_str += export_mesh_edp()

    matrix_types = [ stiffness, Grammian, boundary_Grammian(1,2) ]

    for matrix_type in matrix_types :
        edp_str += export_matrix_edp( **matrix_type )

    FreeFem_output = run_FreeFem( edp_str )

    # print(FreeFem_output)

    mesh = FreeFem_str_to_mesh( FreeFem_output, simple_boundaries = True )

    matrices = {}

    for matrix_type in matrix_types :
        matrices[ matrix_type['matrix_name'] ] = FreeFem_str_to_matrix( FreeFem_output, matrix_type['matrix_name'] )

    for matrix_type in matrix_types :
        print( matrix_type['matrix_name'] )

    epsilon = 1e-4
    M = - matrices[ stiffness['matrix_name'] ] + 1./epsilon*matrices[ boundary_Grammian(1,2)['matrix_name'] ]
    Source = matrices[ Grammian['matrix_name'] ]*np.array( [1]*len( mesh.x ) )
    u = spsolve( M, Source )

    tricontourf( mesh, u )
    # mesh_color =  the_mesh.plot_triangles( triangle_labels = False, alpha = .5 )[0].get_color()
    mesh.plot_boundaries( lw = 2, labels = False )
    # the_mesh.plot_nodes( )

    axis('equal')
    axis('off')

    show()
