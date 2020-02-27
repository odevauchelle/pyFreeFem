# Mesh refinement

## Initial mesh

We first create an unrefined square mesh with FreeFem++.

```python
import pyFreeFem as pyff

script = pyff.edpScript('mesh Th = square(5, 5);')
script += pyff.edpOutput( data_type = 'mesh', name = 'Th' )
Th = script.get_output()['Th']
```
This mesh looks like this:

![Initial mesh](../figures/adaptmesh_0.svg)

## Refinement

We now want to refine the above mesh to accurately represent a two-dimenionnal Gaussian centered at (0.5,0.5). To to so, we use the FreeFem++ function [adaptmesh](https://doc.freefem.org/documentation/mesh-generation.html#the-command-adaptmesh).

```python
# create refinement script

script = pyff.edpScript()
script += pyff.edpInput( name = 'Th', data_type = 'mesh' )
script += 'fespace Vh( Th, P1 );'
script += pyff.edpInput( name = 'u', data_type = 'vector' )
script += 'Th = adaptmesh( Th, u, iso = 1 );'
script += pyff.edpOutput( name = 'Th', data_type = 'mesh' )

# refine Th

for _ in range(3) :
    u = exp( - ( ( Th.x - .5 )**2 + ( Th.y - .5 )**2 )/.1**2 )
    Th = script.get_output( Th = Th, u = u )['Th']

```
The refined mesh looks like this:

![Initial mesh](../figures/adaptmesh_1.svg)
