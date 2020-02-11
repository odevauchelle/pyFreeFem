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
    '''

    def __init__( self, segments, label = None ) :

        self.segments = segments
        self.label = label

    def plot( self, mesh, ax = None, labels = True, **kwargs ) :

        if ax is None :
            ax = gca()

        for segment in self.segments :

            for piece in segment :

                x = mesh.x[ piece ]
                y = mesh.y[ piece ]

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

if __name__ == '__main__' :

    from pylab import *

    x = rand(15)
    y = rand( len(x) )

    the_mesh = TriMesh( x, y )

    the_mesh.plot_triangles( labels = True )
    the_mesh.plot_nodes()

    show()
