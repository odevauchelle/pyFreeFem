# Export mesh to Json

We first create a simple mesh:

```python
import pyFreeFem as pyff

script = pyff.edpScript('''
border Circle( t = 0, 2*pi ){ x = cos(t); y = sin(t); }
mesh Th = buildmesh( Circle(5) );
''')

script += pyff.OutputScript( Th = 'mesh' )

Th = script.get_output()['Th']
```

We then convert the mesh to a Json string:

```python
json_mesh = Th.to_json()

print(json_mesh)
```

```console
>>> {"x": [1.0, 0.309017, 0.309017, -0.809017, -0.809017], "y": [0.0, 0.951057, -0.951057, 0.587785, -0.587785], "triangles": [[0, 1, 2], [1, 3, 2], [3, 4, 2]], "node_labels": [1, 1, 1, 1, 1], "triangle_labels": [0, 0, 0], "boundary_edges": [[0, 0, 1], [1, 0, 1], [2, 0, 1], [2, 1, 1], [0, 2, 1]]}
```

We can then load that string as a dictionary, and create a new mesh from it:

```python
from json import loads

Th = pyff.TriMesh( **loads( json_mesh ) )
```
