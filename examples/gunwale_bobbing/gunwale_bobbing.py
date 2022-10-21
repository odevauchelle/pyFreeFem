#
# O. Devauchelle & Graham Benham, October 2022
#
#
# Benham, G. P., Devauchelle, O., Morris, S. W., & Neufeld, J. A. (2022). Gunwale bobbing. Physical Review Fluids, 7(7), 074804.
#


from pylab import *
from scipy.sparse.linalg import spsolve

import sys
sys.path.append('/home/olivier/git/pyFreeFem/')

import pyFreeFem as pyff

# Mesh


L_box = 10.
H_box = 1.5
L_boat = 1.

mesh_points = dict( boat = [ [L_boat/2, 0] , [-L_boat/2,0] ] )
mesh_points['left_surface'] = [ mesh_points['boat'][-1], [-L_box/2,0] ]
mesh_points['left_wall'] = [ mesh_points['left_surface'][-1], [-L_box/2,-H_box] ]
mesh_points['bottom'] = [ mesh_points['left_wall'][-1], [L_box/2,-H_box] ]
mesh_points['right_wall'] = [ mesh_points['bottom'][-1], [L_box/2,0] ]
mesh_points['right_surface'] = [ mesh_points['right_wall'][-1], [L_boat/2,0] ]

Th = pyff.TriMesh_from_boundaries( boundaries = list( mesh_points.values() ), labels = list( mesh_points.keys() ) )

for _ in range(3) :
    Th = pyff.adaptmesh( Th, iso = 1, hmax = 0.1 )

Th.plot_triangles( color = 'grey', alpha = .5, lw = .5 )
Th.plot_boundaries( labels = ['boat'], color = 'red' )

axis('scaled')
axis('off')

# legend()

script = pyff.InputScript( Th = Th )

# Finite element space and matrices

boundary_name2int, boundary_int2name = Th.get_boundary_label_conversion()

free_surface_int_labels = [ str( boundary_name2int[ name ] ) for name in ( 'left_surface', 'right_surface' ) ]
walls_int_labels = [ str( boundary_name2int[ name ] ) for name in ( 'left_wall', 'right_wall' ) ]

script += pyff.edpScript('fespace Vh( Th, P1 );')
script += pyff.VarfScript(
    stiffness = 'int2d(Th)( dx(u)*dx(v) +  dy(u)*dy(v) )',
    Gramian = 'int2d(Th)( u*v )',
    FreeSurfaceGramian = 'int1d(Th,' + ','.join( free_surface_int_labels ) +  ' )( u*v )',
    WallsGramian = 'int1d(Th,' + ','.join( walls_int_labels ) +  ' )( u*v )',
    BoatGramian = 'int1d(Th,' + str( boundary_name2int['boat'] ) + ' )( u*v )',
    )

ff_out = script.get_output()

k = 3
omega = sqrt( k*tanh( k*H_box ) ) #dispersion relation

M = ff_out['stiffness'] - omega**2*ff_out['FreeSurfaceGramian'] + 1j*k*ff_out['WallsGramian']
B = 1j*omega*ff_out['BoatGramian']*( Th.x + 0 )

phi = spsolve( M, B )

tricontourf( Th, real( phi ) )
show()
