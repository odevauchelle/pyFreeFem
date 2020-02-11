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
///////////////////////////////////////////////////////
'''

###################################
#
# EXPORT MESH
#
###################################

def add_flags( edp_str, flags ) :

    if type( flags ) == type( '' ) :
        flags = 2*[flags]

    start_flag, end_flag = flags

    return 'cout << "# ' + start_flag + '" << endl;\n' + edp_str + 'cout << "# ' + end_flag + '" << endl;\n'



export_nodes =  add_flags( '''
for (int nv = 0; nv < Th.nv; nv++ ) // iteration over all nodes
	{
	cout << Th(nv).x << " " << Th(nv).y << " " << Th(nv).label << endl;
	}
''', 'NODES' )

export_triangles =  add_flags( '''
for (int nt = 0; nt < Th.nt; nt++ ) // iteration over all triangles
	{
	cout << Th[nt][0] << " " << Th[nt][1] << " " << Th[nt][2] << " " << Th[nt].label << endl;
	}
''', 'TRIANGLES' )

export_boundaries =  add_flags( '''
for (int ne = 0; ne < Th.nbe; ne++ ) // iteration over all boundary segments
	{
	cout << Th.be(ne)[0] << " " << Th.be(ne)[1] << " " << Th.be(ne).label << endl;
	}
''', 'BOUNDARIES' )

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

def export_matrix_edp( **kwargs ) :

    edp_str = '''
    varf Vmatrix_name( func, test_func ) = variational_formulation ;
    matrix Mmatrix_name = Vmatrix_name( Vh, Vh ) ;
    Mmatrix_name.resize( func[].n, func[].n );
    '''
    edp_str += add_flags('''
    cout << Mmatrix_name;
    ''', 'MATRIX ' + 'matrix_name' )

    for key in kwargs.keys() :
        edp_str = edp_str.replace( key, kwargs[key] )

    return edp_str.join( [separator]*2 )

###################################
#
# TYPICAL MATRICES
#
###################################

stiffness = dict(
    matrix_name = 'stiffness',
    func = 'u',
    test_func = 'v',
    variational_formulation = 'int2d( Th )( dx(u)*dx(v) + dy(u)*dy(v) )'
    )

Grammian = dict(
    matrix_name = 'Grammian',
    func = 'u',
    test_func = 'v',
    variational_formulation = 'int2d( Th )( u*v )'
    )

def boundary_Neumann( *boundary_label ) :

    boundary_label = list( map( lambda x : str(x), boundary_label ) )

    return dict(
        matrix_name = 'NeumannB' + 'B'.join(boundary_label),
        func = 'u',
        test_func = 'v',
        variational_formulation = 'int1d( Th, ' + ', '.join(boundary_label) + ' )( N.x*dx(u)*v + N.y*dy(u)*v )'
        )

def boundary_Grammian( *boundary_label ) :

    boundary_label = list( map( lambda x : str(x), boundary_label ) )

    return dict(
        matrix_name = 'GrammianB' + 'B'.join(boundary_label),
        func = 'u',
        test_func = 'v',
        variational_formulation = 'int1d( Th, ' + ', '.join(boundary_label) + ' )( u*v + u*v )'
        )



if __name__ == '__main__' :
    print( export_matrix_edp( **stiffness ) )
    print( export_matrix_edp( **Grammian ) )
    print( export_matrix_edp( **boundary_Neumann( 1, 3 ) ) )
    print( export_matrix_edp( **boundary_Grammian( 1, 'circle' ) ) )
