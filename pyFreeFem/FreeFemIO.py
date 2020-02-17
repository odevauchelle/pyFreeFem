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
from tempfile import NamedTemporaryFile

from .TriMesh import TriMesh #, triangle_edge_to_node_edge
from .meshTools.segments import triangle_edge_to_node_edge
from .FreeFemTools.FreeFemStatics import *

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

    boundary_edges = {}

    for edge in FreeFem_mesh['boundaries'] :
        boundary_edges.update( FreeFem_edge_to_boundary_edge( edge, triangles ) )

    mesh = TriMesh(
        x, y,
        triangles = triangles,
        node_labels = node_labels,
        triangle_labels = triangle_labels,
        boundary_edges = boundary_edges
        )

    return mesh



# def FreeFem_to_boundaries( seg_list, reorder = True ) :
#     '''
#     Turns a list of labelled segments into a list of boundaries.
#     '''
#
#     boundaries = []
#
#     for label in set( seg_list[:,2] ) :
#
#         segments = seg_list[ np.where( seg_list[:,2] == label ), :2 ]
#
#         if reorder :
#             segments = reorder_boundary( segments )
#
#         boundaries += [ Boundary( segments = segments, label = label ) ]
#
#     return boundaries

def FreeFem_edge_to_boundary_edge( FreeFem_edge, triangles ) :
    '''
    ( start_node, end_node, label_integer ) -> { ( triangle_index, triangle_node_index ) : label_integer }
    '''

    FreeFem_edge = [ FreeFem_edge[0], FreeFem_edge[1], FreeFem_edge[2] ]

    triangle_index = None

    for nodes_order in [ [0,1], [1,2], [-1,0] ] :
        try :
            triangle_index = triangles[ :, nodes_order].tolist().index( FreeFem_edge[:2] )
            break
        except :
            pass

    if triangle_index is None :

        print('Could not find boundary edge ' + str(FreeFem_edge) + ' in mesh.')
        return {}

    else :

        node_index_in_triangle = triangles[triangle_index].tolist().index( FreeFem_edge[0] )
        return { ( triangle_index, node_index_in_triangle ) : FreeFem_edge[-1]  }
#
# def boundary_edges_to_FreeFem_edges( boundary_edges,  ) :
#     '''
#     { ( triangle_index, triangle_node_index ) : label } -> ( start_node, end_node, label_integer )
#     '''
#
#     if label_to_int is None :
#         label_to_int_func = lambda label : int(label)
#
#     elif type(label_to_int) is type( {} ) :
#         label_to_int_func = lambda label : label_to_int[ label ]
#
#     else :
#         label_to_int_func = label_to_int
#
#     FreeFem_edges = []
#
#     for key in boundary_edges.keys() :
#
#         start_node_index, end_node_index = triangle_edge_to_node_edge( key, triangles )
#
#         FreeFem_edges += [ [ start_node_index + 1, end_node_index + 1, label_to_int_func( boundary_edges[key] ) ] ]
#
#     return np.array( FreeFem_edges )

def savemesh( mesh, filename ) :
    '''
    Saves mesh in FreeFem++ format in a .msh file.
    '''

    with open( filename, 'w' ) as the_file :

        # nv, nt, ne
        the_file.write( str( [ len( mesh.x), len( mesh.triangles ), len( mesh.boundary_edges ) ] )[1:-1].replace(',',' ') + '\n' )

        # vertices
        for node_index in range( len( mesh.x ) ) :
            the_file.write( str( mesh.x[node_index] ) + ' ' + str( mesh.y[node_index] ) + ' ' + str( mesh.node_labels[node_index] ) + '\n' )

        # triangles
        for tri_index, triangle in enumerate( mesh.triangles ) :
            the_file.write( str( list( np.array( triangle ) + 1 ) )[1:-1].replace(',',' ') + ' ' + str( mesh.triangle_labels[tri_index] ) + '\n' )

        # edges
        for edge in mesh.get_boundary_edges() :
            the_file.write( str( list( np.array( edge ) + np.array([ 1, 1, 0 ]) ) )[1:-1].replace(',',' ') + '\n' )

    return filename


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

    with subprocess.Popen( [command], stdout = subprocess.PIPE, stderr = subprocess.PIPE, shell = True, executable="/bin/bash" ) as proc :

        output, error = proc.communicate()

        if not error is b'':
            print('FreeFem++ Error')
            print('---------------')
            print(error)
            print('---------------')

        return output.decode('utf-8')


if __name__ == '__main__' :

    from pylab import *
    from scipy.sparse.linalg import spsolve

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
