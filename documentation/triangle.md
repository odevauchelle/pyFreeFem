# Constrained Delaunay triangulation with the [`triangle`](https://rufat.be/triangle/index.html) library

The [`triangle`](https://rufat.be/triangle/index.html) library can handle constrained Delaunay triangulations. The triangulation it produces can be imported into pyFreeFem.

## Build a [`triangle`](https://rufat.be/triangle/index.html) triangulation

Let's triangulate a box with a river cutting through it:

```python
import triangle as tr
from pylab import *

Ts = dict(
    vertices = [ [-1,-1], [1,-1], [1,1], [-1,1], [0.6,-1], [0.8,0] ],
    segments = [ [ 0, 4 ], [ 4, 1 ], [1, 2], [2,3], [3, 0], [4,5] ]
    )

T = tr.triangulate( Ts, 'pa' )

tr.compare( plt, Ts, T )
```

The result looks like this:

![Compare](./../figures/compare.svg)

## Convert the triangulation into a `TriMesh` object

A [`triangle`](https://rufat.be/triangle/index.html) triangulation is easily imported into pyFreeFem:

```python
import pyFreeFem as pyff

Th = pyff.triangle_to_TriMesh( T )
```

However, the [`triangle`](https://rufat.be/triangle/index.html) does not keep track of the boundary handles. We thus need to add them by hand:

```python
Th.add_boundary_edges( Ts['segments'][:-1], 'basin' )
Th.add_boundary_edges( Ts['segments'][-1], 'river' )
```

We may now refine the mesh and plot it:

```python
for _ in range(3):
    Th = pyff.adaptmesh( Th, hmax = .2, iso = 1 )

Th.plot_triangles()
Th.plot_boundaries()
legend()
```

![Mesh with river](./../figures/river.svg)