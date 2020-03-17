# Input

## Direct input

Let us enter a real number as a FreeFem++ variable:
```python
import pyFreeFem as pyff

script = pyff.InputScript( b = 5.8 )
print( script.get_edp() )
```
The resulting FreeFem++ script looks like
```c++
real b;
cin >> b;
```
By default, this declares the variable `b`, a uses the `cin` input of FreeFem++ to populate it.

## Postponed input

The variable's value needs not be specified, but then its type is required:
```python
import pyFreeFem as pyff

script = pyff.InputScript( b = 'real' )
script += 'cout << b << endl;'
print( script.run( b = 6.7 ) )
print( script.run( b = 3 ) )
```

## Multiple inputs

Inputs can be of different types, postponed or not:
```python
import pyFreeFem as pyff

script = pyff.InputScript( a = 5, b = 'real' )
script += 'cout << a*b << endl;'
print( script.run( b = 6.7 ) )
```

## Declaration

FreeFem++ can't handle multiple declarations, so the following fails:
```python
import pyFreeFem as pyff

script = pyff.InputScript( a = 5 )
script += pyff.InputScript( a = 14 )
script.run()
```
The solution is to drop the variable's declaration:
```python
script += pyff.InputScript( a = 14, declare = False )
```
This, however, affects the entire call of `InputScript`, that is, none of its arguments will be declared.

# Output

Ouputs are somewhat simpler than inputs, but the output variable's type needs to be declared:
```python
script = pyff.InputScript( a = 5.1 )
script += pyff.OutputScript( a = 'real' )
print( script.get_output()['a'] )
```
The function `get_output` parses the output of FreeFem++, and returns the result as a dictionnary.

# I/O Types

pyFreeFem currently handles the follwing variable types:

| FreeFem++      | pyFreeFem  |  Python |
| :------------: | :----------: | :----------: |
| `int` | `int`  | `int`  |
| `real` | `real`  | `float`  |
| `mesh` | `mesh` | `TriMesh` |
| `matrix` | `matrix` | `csr_matrix` |
| fespace name (typically `Vh`) | `vector` | `list`, `array`, `ndarray` |
