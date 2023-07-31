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


from .edpTools import flagize

separator = '''
///////////////////////////////////////////////////////\n
'''

default_variable_names = {
    '_u_' : 'u',
    '_v_' : 'v',
    '_VhU_' : 'Vh',
    '_VhV_' : 'Vh',
    '_Th_' : 'Th'
    }



###################################
#
# EXPORT MESH to FreeFem++
#
###################################

def add_flags( edp_str, flags ) :

    if type( flags ) == type( '' ) :
        flags = 2*[flags]

    start_flag, end_flag = flags

    edp_str = 'cout << "' + start_flag + '" << endl;\n' + edp_str + 'cout << "' + end_flag + '" << endl;\n'

    return edp_str

export_nodes =  add_flags( '''
for (int nv = 0; nv < _Th_.nv; nv++ )
	{
	cout << _Th_(nv).x << " " << _Th_(nv).y << " " << _Th_(nv).label << endl;
	}
''', flagize( 'nodes' ) )

export_triangles =  add_flags( '''
for (int nt = 0; nt < _Th_.nt; nt++ )
	{
	cout << _Th_[nt][0] << " " << _Th_[nt][1] << " " << _Th_[nt][2] << " " << _Th_[nt].label << endl;
	}
''', flagize( 'triangles' ) )

export_boundaries =  add_flags( '''
for (int ne = 0; ne < _Th_.nbe; ne++ )
	{
	cout << _Th_.be(ne)[0] << " " << _Th_.be(ne)[1] << " " << _Th_.be(ne).label << endl;
	}
''', flagize( 'boundaries' ) )


def export_mesh_edp( **kwargs ) :

    edp_str = separator.join( [ export_nodes, export_triangles, export_boundaries ] )

    for name in kwargs.items() :
        edp_str = edp_str.replace( *name )

    return edp_str

###################################
#
# EXPORT MATRICES
#
###################################

# def create_varf_matrix( **kwargs ) :
#
#     edp_str = '''
#     varf V_matrix_name_( _u_, _v_ ) = _variational_expression_ ;
#     matrix M_matrix_name_ = V_matrix_name_( _Vh_u_, _Vh_v_ ) ;
#     '''
#
#     for name in kwargs.items() :
#         edp_str = edp_str.replace( *name )
#
#     return  edp_str

def export_matrix_edp( create_and_add_flags = True, **kwargs ) :

    edp_str = ''
    # M_matrix_name_.resize( _u_[].n, _v_[].n );

    export_edp_str = '''
    cout << M_matrix_name_;
    '''

    if create_and_add_flags :
        try :
            flag = flagize( kwargs['_matrix_name_'] )
        except :
            flag = flagize( '_matrix_name_' )

        export_edp_str += add_flags( export_edp_str, flag )

    edp_str += export_edp_str

    for name in kwargs.items() :
        edp_str = edp_str.replace( *name )

    return edp_str

###################################
#
# EXPORT VECTORS
#
###################################

def export_vector_edp( **kwargs ) :

    # edp_str = 'cout << base_func[];\n'

    edp_str = '''
    for (int nVector = 0; nVector < _u_.n; nVector++ )
    	{
    	cout << _u_[][nVector] << endl;
    	}
    '''

    for name in kwargs.items() :
        edp_str = edp_str.replace( *name )

    return edp_str
