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
from pylab import gca, mean, array

TriMesh_structure ='''
TriMesh attributes:

    - Nodes coordinates (from Triangulation)
    - Oriented triangles (from Triangulation)
    - Nodes labels (list)
    - Triangle labels (list)

    - Oriented edges with labels (dict)

        - Triangle index
        - Node index (0,1 or 2)
        - Label

TriMesh inner classes:

    - Boundary

    Boundary attributes:

        - Label (that of edges)
        - List of segments
            - oriented list of nodes (connected by edges)
'''

class TriMesh( mptri.Triangulation ) :

    '''
    Triangular mesh structure suited for FreeFem++.
    This class is inherited from matplotlib.tri.Triangulation.
    Additionnal features are:
        Node labels;
        Import and export from and to FreeFem mesh.
    '''

    def __init__( self, x, y, triangles = None, boundary_edges = None, node_labels = None, triangle_labels = None ) :

        '''
        Create a triangular mesh suited for FreeFem++.

        Parameters :
        ----------

        x, y : n*1 lists of floats
            Nodes coordinates
        triangles : 3*n list of integers, optional
            Triangles indices
        '''

        super( TriMesh, self ).__init__( x, y, triangles )


        if node_labels is None :
            node_labels = [0]*len( self.x )

        self.node_labels = node_labels

        if triangle_labels is None :
            triangle_labels = [0]*len( self.triangles )

        self.triangle_labels = triangle_labels

        if boundary_edges is None :
            boundary_edges = {}

        self.boundary_edges = boundary_edges

    def plot_triangles( self, labels = None, ax = None, **kwargs ) :

        if ax is None :
            ax = gca()

        triangles_plot = ax.triplot( self, **kwargs )

        label_style = dict(va = 'center', ha = 'center', color =  triangles_plot[0].get_color())

        if not labels is None :

            for i, triangle in enumerate( self.triangles ) :

                x, y = mean( self.x[triangle] ), mean( self.y[triangle] )

                if labels is 'label' :
                    ax.text( x, y , self.triangle_labels[i], **label_style )

                elif labels is 'index' :
                    ax.text( x, y , i, **label_style )

        return triangles_plot

    def plot_nodes( self, labels = None, ax = None, **kwargs ) :

        if ax is None :
            ax = gca()

        nodes_plot_style = dict( marker = 'o', linestyle = 'None', ms = 12 )
        nodes_plot_style.update( kwargs )

        nodes_plot = ax.plot( self.x, self.y, **nodes_plot_style )

        if not labels is None :

            label_style = dict(va = 'center', ha = 'center', color = 'w')

            for i in range( len( self.x ) ) :

                if labels is 'label' :
                    ax.text( self.x[i], self.y[i], self.node_labels[i], **label_style )

                elif labels is 'index' :
                    ax.text( self.x[i], self.y[i], i, **label_style )

        return nodes_plot

    def plot_edges( self, labels = None, ax = None, **kwargs ) :

        if ax is None :
            ax = gca()

        for key in self.edges.keys() :

            nodes = 

            nodes_plot = ax.plot( self.x, self.y, **nodes_plot_style )


if __name__ == '__main__' :

    from pylab import *

    x = rand(15)
    y = rand( len(x) )

    boundary_edges = { ( 0, 0 ) : 'B'  }

    mesh = TriMesh( x, y, boundary_edges = boundary_edges )


    color = mesh.plot_triangles( labels = 'index' )[0].get_color()
    mesh.plot_nodes( labels = 'index', color = color )

    show()
