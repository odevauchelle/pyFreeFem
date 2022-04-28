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
from pylab import gca, mean, array, nan

if __name__ == '__main__' :
    from meshTools.segments import *

else :
    from .meshTools.segments import *

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
        self.boundary_edges = {}

        if not boundary_edges is None :
            self.add_boundary_edges( boundary_edges )


    def add_boundary_edges( self, boundary_edges, label = None ) :
        '''
        Add a boundary.
        '''

        try :
            # assume { (triangle_index, node_in_triangle) : label, ...  }
            boundary_edges.keys() # check if dictionnary before updating
            self.boundary_edges.update( boundary_edges )

        except :

            try :
                # assume [ [ triangle_index, node_in_triangle, label ], ... ]
                boundary_edges[0][2] # check format; label value implies triangle-based indexing
                self.add_boundary_edges( edges_to_boundary_edges( boundary_edges ) )

            except :

                if label is None :
                    label = invent_label( self.boundary_edges.values() )

                try :
                    # assume [ [i,j], [k,l], ... ] where i, j, k, l are node node_indices
                    boundary_edges[0][1] # check format;
                    self.add_boundary_edges( node_index_to_triangle_index_edges( boundary_edges, self.triangles, label = label ) )

                except :
                    # assume [ first_node, second_node, ... ]
                    edges = [ [ boundary_edges[i], boundary_edges[i + 1] ] for i in range( len(boundary_edges) - 1 ) ]

                    self.add_boundary_edges( edges, label = label )


    def get_boundary_label_conversion( self ) :
        return label_conversion( list( self.boundary_edges.values() ) )

    def rename_boundary( self, new_names, verbose = False ) :

        for edge in self.boundary_edges.keys() :
            try :
                self.boundary_edges[edge] = new_names[ self.boundary_edges[edge] ]
            except :
                if verbose :
                    print( 'No new name for ' + str( self.boundary_edges[edge] ) )
                pass

    def get_boundary_edges( self, label_type = 'int', index_type = 'node' ) :

        '''
        Method to convert boundary segments into an array.

        Typically:
        { ( triangle_index, node_index_in_triangle ) : raw_label }
            -> [ start_node_index, end_node_index, int_label ]
        '''

        if label_type == 'int' :
            label_to_int, _ = self.get_boundary_label_conversion()

        edges = []

        for edge in self.boundary_edges.keys() :

            if index_type == 'node' :
                edge_indices = list( triangle_edge_to_node_edge( edge, self.triangles ) )

            elif index_type == 'triangle' :
                edge_indices = list( edge )

            if label_type == 'int' :
                label = label_to_int[ self.boundary_edges[ edge ] ]

            else :
                label = self.boundary_edges[ edge ]

            edges += [ edge_indices + [ label ] ]

        return edges

    def get_boundaries( self ) :

        '''
        A boundary is a dictionnary indexed by label.
        Each value is a list of segments.
        Each segment is an oriented list of nodes
        '''

        boundaries = {}

        for edge in self.get_boundary_edges( label_type = 'raw' ) :

            label = edge[-1]

            try :
                edges = boundaries[label]
            except :
                edges = []

            edges += [ edge[:-1] ] # each edge is properly oriented, but edges are not sorted at this point

            boundaries.update( { label : edges } )

        for label in boundaries.keys() :

            boundaries[label] = edges_to_segments(  boundaries[label] )

        return boundaries


    def plot_triangles( self, labels = None, ax = None, **kwargs ) :

        if ax is None :
            ax = gca()

        triangles_plot = ax.triplot( self, **kwargs )

        label_style = dict(va = 'center', ha = 'center', color =  triangles_plot[0].get_color())

        if not labels is None :

            for i, triangle in enumerate( self.triangles ) :

                x, y = mean( self.x[triangle] ), mean( self.y[triangle] )

                if labels == 'label' :
                    ax.text( x, y , self.triangle_labels[i], **label_style )

                elif labels == 'index' :
                    ax.text( x, y , i, **label_style )

        return triangles_plot

    def plot_nodes( self, labels = None, ax = None, **kwargs ) :

        if ax is None :
            ax = gca()

        nodes_plot_style = dict( marker = 'o', linestyle = 'None', ms = 12 )
        nodes_plot_style.update( kwargs )

        nodes_plot = ax.plot( self.x, self.y, **nodes_plot_style )

        if not labels is None :

            label_style = dict( va = 'center', ha = 'center', color = 'w' )

            for i in range( len( self.x ) ) :

                if labels == 'label' :
                    ax.text( self.x[i], self.y[i], self.node_labels[i], **label_style )

                elif labels == 'index' :
                    ax.text( self.x[i], self.y[i], i, **label_style )

        return nodes_plot

    def plot_boundaries( self, labels = None, ax = None, **kwargs ) :

        if ax is None :
            ax = gca()

        if labels is None :
            labels = self.get_boundaries().keys()

        try:
            enumerate(labels)
        except:
            labels = [labels] # assume single label

        for label in labels :

            segments = self.get_boundaries()[label]

            boundary_kwargs = { 'label' : label }
            boundary_kwargs.update( kwargs )

            x = []
            y = []

            for segment in segments :

                x += [nan] + list( self.x[segment] )
                y += [nan] + list( self.y[segment] )

            ax.plot(x, y, **boundary_kwargs )

    def plot_edges( self, labels = None, ax = None, **kwargs ) :

        if ax is None :
            ax = gca()

        for edge in self.boundary_edges.keys() :

            node_indices = list( triangle_edge_to_node_edge( edge, self.triangles ) )

            x, y = self.x[node_indices], self.y[node_indices]

            edge_plot = ax.plot( x, y, **kwargs )

            if labels == 'label' :
                label_style = dict(va = 'center', ha = 'center', color = edge_plot[0].get_color() )
                ax.plot( mean(x), mean(y), 'ow', ms = 12 )
                ax.text( mean(x), mean(y), self.boundary_edges[edge], **label_style )



if __name__ == '__main__' :

    from pylab import *

    x = array([1,5,6,3,4])
    y = array([1,4,7,3,0])

    boundary_edges = { ( 3, 0 ) : 'B',  ( 3, 1 ) : 'B' }

    mesh = TriMesh( x, y, boundary_edges = boundary_edges )
    color = mesh.plot_triangles( labels = 'index' )[0].get_color()
    mesh.plot_nodes( labels = 'index', color = color )
    mesh.plot_edges( labels = 'label', color = 'red' )

    print( mesh.get_boundaries() )

    show()
