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

separator = '''
///////////////////////////////////////////////////////\n
'''

from .edpTools import flagize

###################################
#
# EXPORT MESH
#
###################################

def add_flags( edp_str, flags ) :

    if type( flags ) == type( '' ) :
        flags = 2*[flags]

    start_flag, end_flag = flags

    edp_str = 'cout << "' + start_flag + '" << endl;\n' + edp_str + 'cout << "' + end_flag + '" << endl;\n'

    return edp_str

export_nodes =  add_flags( '''
for (int nv = 0; nv < Th.nv; nv++ )
	{
	cout << Th(nv).x << " " << Th(nv).y << " " << Th(nv).label << endl;
	}
''', flagize( 'nodes' ) )

export_triangles =  add_flags( '''
for (int nt = 0; nt < Th.nt; nt++ )
	{
	cout << Th[nt][0] << " " << Th[nt][1] << " " << Th[nt][2] << " " << Th[nt].label << endl;
	}
''', flagize( 'triangles' ) )

export_boundaries =  add_flags( '''
for (int ne = 0; ne < Th.nbe; ne++ )
	{
	cout << Th.be(ne)[0] << " " << Th.be(ne)[1] << " " << Th.be(ne).label << endl;
	}
''', flagize( 'boundaries' ) )

# after: https://www.ljll.math.upmc.fr/pipermail/freefempp/2008/000164.html

'''
int[int] AllLabels = labels(Th);
'''


def export_mesh_edp( **kwargs ) :

    edp_str = separator.join( [ export_nodes, export_triangles, export_boundaries ] )

    for key in kwargs.keys() :
        edp_str = edp_str.replace( key, kwargs[key] )

    return edp_str

###################################
#
# EXPORT MATRICES
#
###################################

def create_varf_matrix( **kwargs ) :

    edp_str = '''
    varf Vmatrix_name( base_func, test_func ) = variational_formulation ;
    matrix Mmatrix_name = Vmatrix_name( Vh, Vh ) ;
    '''

    for key in kwargs.keys() :
        edp_str = edp_str.replace( key, kwargs[key] )

    return  edp_str

def export_matrix_edp( create_varf = True, create_and_add_flags = True, **kwargs ) :

    edp_str = ''

    if create_varf :
        edp_str += create_varf_matrix()
        edp_str += separator

    export_edp_str = '''
    Mmatrix_name.resize( base_func[].n, base_func[].n );
    cout << Mmatrix_name;
    '''

    if create_and_add_flags :
        try :
            flag = flagize( kwargs['matrix_name'] )
        except :
            flag = flagize( 'matrix_name' )

        export_edp_str += add_flags( export_edp_str, flag )

    edp_str += export_edp_str

    for key in kwargs.keys() :
        edp_str = edp_str.replace( key, kwargs[key] )

    return edp_str

###################################
#
# EXPORT VECTORS
#
###################################

def export_vector_edp( **kwargs ) :

    # edp_str = 'cout << base_func[];\n'

    edp_str = '''
    for (int nVector = 0; nVector < base_func.n; nVector++ )
    	{
    	cout << base_func[][nVector] << endl;
    	}
    '''


    for key in kwargs.keys() :
        edp_str = edp_str.replace( key, kwargs[key] )

    return edp_str

###################################
#
# TYPICAL MATRICES
#
###################################

stiffness = dict(
    matrix_name = 'stiffness',
    base_func = 'u',
    test_func = 'v',
    variational_formulation = 'int2d( Th )( dx(u)*dx(v) + dy(u)*dy(v) )'
    )

Grammian = dict(
    matrix_name = 'Grammian',
    base_func = 'u',
    test_func = 'v',
    variational_formulation = 'int2d( Th )( u*v )'
    )

def boundary_Neumann( *boundary_label ) :

    boundary_label = list( map( lambda x : str(x), boundary_label ) )

    return dict(
        matrix_name = 'NeumannB' + 'B'.join(boundary_label),
        base_func = 'u',
        test_func = 'v',
        variational_formulation = 'int1d( Th, ' + ', '.join(boundary_label) + ' )( N.x*dx(u)*v + N.y*dy(u)*v )'
        )

def boundary_Grammian( *boundary_label ) :

    boundary_label = list( map( lambda x : str(x), boundary_label ) )

    return dict(
        matrix_name = 'GrammianB' + 'B'.join(boundary_label),
        base_func = 'u',
        test_func = 'v',
        variational_formulation = 'int1d( Th, ' + ', '.join(boundary_label) + ' )( u*v + u*v )'
        )



if __name__ == '__main__' :
    print( export_mesh_edp() )
    # print( export_matrix_edp( **stiffness ) )
    # print( export_matrix_edp( **Grammian ) )
    # print( export_matrix_edp( **boundary_Neumann( 1, 3 ) ) )
    # print( export_matrix_edp( **boundary_Grammian( 1, 'circle' ) ) )
