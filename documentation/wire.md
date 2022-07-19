# Electric field around a resistive wire

We want to calculate the electric field around an electric wire. For this, we need to solve the Laplace equation over a domain with a one-dimensional boundary cutting into it.

## Domain and mesh

We first define a square box and a wire:

```python
from pylab import *
import pyFreeFem as pyff

box_points = [ [ .5, 0 ], [.5,1], [-.5,1], [-.5,0] ]
wire_points = [ [0,0],[0,.2],[.1,.4] ]
```

We then build the mesh, and define the two boundaries.

```python
Th = pyff.TriMesh( *array( box_points + wire_points ).T )
Th.add_boundary_edges( range( len( box_points ) ) , 'box' )
Th.add_boundary_edges( range( len( box_points ), len( box_points ) + len( wire_points ) ) , 'wire' )
```

Here is the result:

```python
Th.plot_triangles( color = 'grey', labels = 'index' )
Th.plot_boundaries()
Th.plot_nodes( color = 'grey', labels = 'index' )

legend(loc = 'upper center')

axis('scaled')
axis('off')
show()
```

!['Mesh around wire'](../figures/wire_mesh.svg)

## Matrices

To create the finite-elements matrices, we first need to import the mesh in a pyFreeFem script, and create the associated P1 element space:

```python
script = pyff.InputScript( Th = Th )
script += pyff.edpScript('fespace Vh( Th, P1 );')
```

We then write the variational forms we're interested in. First, the stiffness matrix:

```python
variational_forms = dict( stiffness = 'int2d(Th)( dx(u)*dx(v) +  dy(u)*dy(v) )' )
```

We now want to calculate the Gramian on each boundary. Let's first convert the boundary names into FreeFem indexes:

```python
name_to_index, _ = Th.get_boundary_label_conversion()
```

We can then use the resulting dictionary to create the Gramians:

```python
for name, index in name_to_index.items() :
    variational_forms.update( { 'BoundaryGramian_' + name : 'int1d(Th,' + str(index) + ')( u*v )' } )
```

```console
>>> {'stiffness': 'int2d(Th)( dx(u)*dx(v) +  dy(u)*dy(v) )', 'BoundaryGramian_box': 'int1d(Th,1)( u*v )', 'BoundaryGramian_wire': 'int1d(Th,2)( u*v )'}
```

Finally, we compute the matrices:

```python
script += pyff.VarfScript( **variational_forms )
matrices = script.get_output()
```

## A look at the matrices

Here's how our matrices look like:

```python
fig, axs = subplots( ncols = len( matrices ), figsize = ( 8, 3 ), sharey = True )

for i, ( name, mat ) in enumerate( matrices.items() ) :

    ax = axs[i]
    ax.imshow( mat.toarray() )
    ax.set_title(name)
```

!['Matrices'](../figures/wire_matrices.svg)

The third matrix is the Gramian along the wire. As expected, its only non-zero coefficients are those which involve nodes 4, 5 and 6.

## Absorbing boundary condition

As a warm-up, we can assume that the wire is a perfect conductor (absorbing boundary), and solve the Laplace equation around it. Let's call the potential $v$, and assume that it is 1 on the box, and 0 along the wire. These conditions translate into $v-1=\epsilon \partial_n v$ on the box, and $v=\epsilon \partial_n v$ on the wire.

We further assume that the bottom line (which we haven't defined as a boundary) is reflective. Altogether, the weak formulation of our problem reads ($\hat{v}$ is the test function):

$$
\iint \nabla \hat{v} \cdot \nabla v + \dfrac{1}{\epsilon} \oint \hat{v} v - \dfrac{1}{\epsilon} \int_{\mathrm{box}} \hat{v} v
$$

where the reflective boundary does not appear. In matrix form, the above problem reads

$$
\left( \mathbf{S} + \dfrac{1}{\epsilon} \left( \mathbf{G}_{\mathrm{box}} + \mathbf{G}_{\mathrm{wire}} \right) \right) \cdot V = \dfrac{1}{\epsilon} \mathbf{G}_{\mathrm{box}} \cdot \mathbb{1}
$$

We can use `spsolve` to solve this linear problem:


```python
from scipy.sparse.linalg import spsolve

epsilon = 1e-6
ones_vector = Th.x*0 + 1.

v = spsolve(
    matrices['stiffness'] + 1/epsilon*( matrices['BoundaryGramian_box'] + matrices['BoundaryGramian_wire'] ),
    matrices['BoundaryGramian_box']*ones_vector
    )
```

Here is how the result looks like

```python
figure()
ax_v = gca()

ax_v.tricontourf( Th, v )

Th.plot_boundaries(ax = ax_v, clip_on = False)
Th.plot_triangles(ax = ax_v, color = 'w', alpha = .3 )
```

!['Absorbing boundary'](../figures/wire_field.svg)
