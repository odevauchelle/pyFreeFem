# PyFreeFem documentation

## Finite element matrices

To get finite element matrices, we first need to create a mesh, and define a finite element space. The VarfBlock function then creates, and exports, the matrix corresponding to a variational formulation:

```python
import pyFreeFem as pyff

script = pyff.edpScript('''
border Circle( t = 0, 2*pi ){ x = cos(t); y = sin(t); }
mesh Th = buildmesh( Circle(20) );
fespace Vh( Th, P1 );
Vh u,v;
''' )

# Create and export stiffness matrix
script += pyff.VarfScript( stiffness = 'int2d(Th)( dx(u)*dx(v) +  dy(u)*dy(v) )')
stiffness = script.get_output()['stiffness']
print(stiffness)
```

```console
>>> (0, 0)	1.72789742056
>>> (0, 1)	-0.422931342447
>>> (0, 8)	-0.331359535237
>>> (0, 9)	-0.973606542874
...
```

## Solve a simple problem

We first create a mesh, and import the associated matrices.

```python
import pyFreeFem as pyff
import matplotlib.pyplot as pp

script = pyff.edpScript('''
real smallRadius = .3;
border outerCircle( t = 0, 2*pi ){ x = cos(t); y = 0.8*sin(t); }
border innerCircle( t = 2*pi, 0 ){ x = .5 + smallRadius*cos(t); y = smallRadius*sin(t); }
mesh Th = buildmesh( outerCircle(100) + innerCircle(40) );

fespace Vh( Th, P1 );
Vh u,v;
''')

script += pyff.OutputScript( Th = 'mesh' )

script += pyff.VarfScript(
    stiffness = 'int2d(Th)( dx(u)*dx(v) +  dy(u)*dy(v) )',
    Gramian = 'int2d(Th)( u*v )',
    boundary_Gramian = 'int1d(Th, 1, 2)( u*v )'
    )

ff_output = script.get_output()
Th = ff_output['Th']

Th.plot_triangles( color = 'k', alpha = .2, lw = .5 )
Th.plot_boundaries()
pp.show()
```
The mesh looks like this:

![Mesh with a hole](../figures/solve.svg)

We now solve Poisson's equation on this mesh, with absorbing boundary conditions on the two boundaries. We do this with spsolve from the scipy library.

```python
from scipy.sparse.linalg import spsolve
import numpy as np

epsilon = 1e-4
M = - ff_output['stiffness'] + 1./epsilon*ff_output['boundary_Gramian']
Source = ff_output['Gramian']*np.array( [1]*len( Th.x ) )
u = spsolve( M, Source )

pp.tricontourf( Th, u )
Th.plot_boundaries( color = 'black' )
pp.show()
```
Here is the result:

![Mesh with a hole](../figures/solve_2.svg)

## Mess with the mesh

We now create a mesh with FreeFem++, import it as a TriMesh, change its boundaries and export it back to FreeFem++. [Exports and imports](./documentation/IO.md) to and from FreeFem++ are what pyFreeFem was written for.
```python
import pyFreeFem as pyff

# Create mesh with FreeFem++
script = pyff.edpScript( '''
    border Circle( t = 0, 2*pi ){ x = cos(t); y = sin(t); }
    mesh Th = buildmesh( Circle(150) );
    ''' )
script += pyff.OutputScript( Th = 'mesh' )
Th = script.get_output()['Th']

# Change mesh
for triangle_index in sample( range( len( Th.triangles ) ), 20 ) :
    Th.boundary_edges.update( { ( triangle_index, 0 ) : 2 } )

Th.rename_boundary( {1:'initial', 2:'new'} )
```
The mesh looks like this:

![Messed up mesh](../figures/mesh_IO_mesh.svg)

We want to solve the Poisson equation on this new mesh. Let us first calculate the finite-elements matrices we need.

```python
# Export mesh back to FreeFem
script = pyff.InputScript( Th = Th )

# calculate FEM matrices
script += '''
    fespace Vh( Th, P1 ) ;
    Vh u, v ;
    '''

matrices = {
    'stiffness' : 'int2d(Th)( dx(u)*dx(v) +  dy(u)*dy(v) )',
    'Grammian' : 'int2d(Th)( u*v )',
    'boundary_Grammian' : 'int1d(Th, 1, 2)( u*v )'
    }

script += pyff.VarfScript( **matrices )
matrices = script.get_output()
```
We may now solve our finite-element problem as [above](#solve-a-simple-problem). Here is the result:

![Poisson on messed up mesh](../figures/mesh_IO_field.svg)

A more useful mesh change, perhaps, is to [refine it](./documentation/adaptmesh.md).

# Detailed examples

- [Input and output](./IO.md) to and from FreeFem++

- [Build](./build_your_own_mesh.md) a mesh from scratch

- [Refine](./adaptmesh.md) a mesh manually

- [Mix](./mixed_FE_spaces.md) finite-element spaces on the same mesh

- Compute [travel times](./travel_time.md) along stream lines
