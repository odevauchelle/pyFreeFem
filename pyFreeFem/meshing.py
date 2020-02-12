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

import matplotlib.tri as mptri
from pylab import gca, mean

class Boundary :
    '''
    A boundary for FreeFem++, that is, a labelled list of segments.
    A segment is a list of egdes.
    An edge is a list ot two vertex indices.
    '''

    def __init__( self, segments, label = None ) :

        self.segments = segments
        self.label = label

    def plot( self, mesh, ax = None, labels = True, **kwargs ) :

        if ax is None :
            ax = gca()

        for segment in self.segments :

            for edge in segment :

                x = mesh.x[ edge ]
                y = mesh.y[ edge ]

                color = ax.plot( x, y, **kwargs )[0].get_color()

                if labels :
                    ax.plot( mean(x), mean(y), 'ow', ms = 11 )#, mec = color )
                    ax.text( mean(x), mean(y), str(self.label), va = 'center', ha = 'center', color = color )

class TriMesh( mptri.Triangulation ) :

    '''
    Triangular mesh structure suited for FreeFem++.
    This class is inherited from matplotlib.tri.Triangulation.
    Additionnal features are:
        Node labels;
        Import and export from and to FreeFem mesh.
    '''

    def __init__( self, x, y, triangles = None, node_labels = None, triangle_labels = None, boundaries = None ) :

        '''
        Create a triangular mesh suited for FreeFem++.

        Parameters :
        ----------

        x, y : n*1 lists of floats
            Nodes coordinates
        triangles : 3*n list of integers, optional
            Triangles indices
        labels : 1*n list ( str, int, etc. ), optional
            Nodes labels. Inner nodes are labelled "0" in FreeFem++. Default is 0.
        boundaries : list of Boundary objects, optional
            A boundary for FreeFem++, that is, a list of node indices.
            This list is ordered, and therefore oriented.
            Each of its its segments should belong to a triangle, although this is not enforced a priori.
        '''

        super( TriMesh, self ).__init__( x, y, triangles )

        if node_labels is None :
            self.node_labels = [0]*len( self.x )
        else :
            self.node_labels = node_labels

        if triangle_labels is None :
            self.triangle_labels = [0]*len( self.triangles )
        else :
            self.triangle_labels = triangle_labels

        if boundaries is None :
            self.boundaries = []
        else :
            self.boundaries = boundaries

    def plot_triangles( self, ax = None, labels = False, **kwargs ) :
        '''
        Plot mesh structure.
        '''

        if ax is None :
            ax = gca()

        triangles_plot =  ax.triplot( self.x, self.y, triangles = self.triangles, **kwargs )

        if labels :

            for i, triangle in enumerate( self.triangles ) :
                x, y = mean( self.x[triangle] ), mean( self.y[triangle] )
                ax.text( x, y , self.triangle_labels[i], va = 'center', ha = 'center', color =  triangles_plot[0].get_color() )

        return triangles_plot

    def plot_nodes( self, ax = None, color = 'tab:blue' ) :
        '''
        Plot mesh nodes with labels.
        '''

        if ax is None :
            ax = gca()

        for i in range( len( self.x ) ) :
            ax.plot( self.x[i], self.y[i], 'o', ms = 10, color = color )
            ax.text( self.x[i], self.y[i], self.node_labels[i], va = 'center', ha = 'center', color = 'w' )

    def plot_boundaries( self, ax = None, labels = False, **kwargs ) :

        '''
        Plot meshboundaries with labels.
        '''

        if ax is None :
            ax = gca()

        for boundary in self.boundaries :

            color = ax.plot( [], [], **kwargs )[0].get_color()

            boundary.plot( self, ax = ax, color = color, labels = labels )

    def get_integer_boundary_labels( self ) :
        '''
        Returns a list of integer labels. Useful for FreeFem++, which cannot load string labels.
        '''

        int_labels = [0]

        for boundary in self.boundaries :

            try :
                int_labels += [ int( boundary.label ) ]
            except :
                int_labels += [ max( int_labels ) + 1 ]

        return int_labels[1:]


    def get_boundary_edges( self, integer_label = True ):
        '''
        Returns a list of labelled egdes. These make up the boundaries.
        An edge is [ index_vertex_1, index_vertex_2, boundary_label ]
        '''

        edges = []

        for nb, boundary in enumerate( self.boundaries ) :

            if integer_label :
                label = self.get_integer_boundary_labels()[nb]
            else :
                label = boundary.label

            for segment in boundary.segments :
                for edge in segment :
                    edges += [ edge + [ label ] ]

        return edges

    def save( self, filename ) :
        '''
        Saves mesh in FreeFem++ format in a .msh file.
        '''

        with open( filename, 'w' ) as the_file :

            # nv, nt, ne
            the_file.write( str( [ len(self.x), len( self.triangles ), len( self.get_boundary_edges() ) ] )[1:-1].replace(',',' ') + '\n' )

            # vertices
            for node_index in range( len( self.x ) ) :
                the_file.write( str( self.x[node_index] ) + ' ' + str( self.y[node_index] ) + ' ' + str( self.node_labels[node_index] ) + '\n' )

            # triangles
            for tri_index, triangle in enumerate( self.triangles ) :
                the_file.write( str( list( array( triangle ) + 1 ) )[1:-1].replace(',',' ') + ' ' + str( self.triangle_labels[tri_index] ) + '\n' )
            # edges
            for edge in self.get_boundary_edges() :
                the_file.write( str( list( array( edge ) + 1 ) )[1:-1].replace(',',' ') + '\n' )

        return filename

if __name__ == '__main__' :

    from pylab import *

    x = rand(15)
    y = rand( len(x) )
    boundaries = [ Boundary( [[[1,2]]], label = 'bound' ), Boundary( [[[4,5]]], label = '3' ) ]

    mesh = TriMesh( x, y, boundaries = boundaries )

    mesh.plot_triangles( labels = True )
    mesh.plot_nodes()
    mesh.plot_boundaries( labels = True, color = 'red' )

    print(mesh.get_boundary_edges())

    print( mesh.save('../sandbox/mesh.msh') )


    show()
