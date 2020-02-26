# Mesh refinement

## Initial mesh

We first create an unrefined square mesh with FreeFem++.

```python
import pyFreeFem as pyff

script = pyff.edpScript( 'mesh Th = square(5, 5);')
script += pyff.edpOutput( type = 'mesh', name = 'Th' )

Th_initial = script.get_output()['Th']
```
This mesh looks like this:

![Initial mesh](../figures/adaptmesh_0.svg)

## Refinement

We now want to refine the above mesh to accurately represent a two-dimenionnal Gaussian centered at (0.5,0.5). To to so, we use the FreeFem++ function [adaptmesh](https://doc.freefem.org/documentation/mesh-generation.html#the-command-adaptmesh).

```python
Th = Th_initial

for _ in range(3) :

    u = exp( - ( ( Th.x - .5 )**2 + ( Th.y - .5 )**2 )/.1**2 )

    script = pyff.edpScript()
    script += pyff.edpInput( Th, name = 'Th' )
    script += 'fespace Vh( Th, P1 );'
    script += pyff.edpInput( u, name = 'u' )
    script +='Th = adaptmesh( Th, u, iso = 1 );'
    script += pyff.edpOutput( name = 'Th', type = 'mesh' )

    Th = script.get_output()['Th']
```
The refined mesh looks like this:

![Initial mesh](../figures/adaptmesh_1.svg)
