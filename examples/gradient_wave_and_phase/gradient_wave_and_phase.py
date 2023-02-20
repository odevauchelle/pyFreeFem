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

#########################
#
# build the mesh
#
#########################

L_box = 6.
H_box = 1.5

mesh_points = {}

mesh_points['surface'] = [ [L_box/2,0], [-L_box/2,0] ]
mesh_points['left_wall'] = [ mesh_points['surface'][-1], [-L_box/2,-H_box] ]
mesh_points['bottom'] = [ mesh_points['left_wall'][-1], [L_box/2,-H_box] ]
mesh_points['right_wall'] = [ mesh_points['bottom'][-1], [L_box/2,0] ]

Th = pyff.TriMesh_from_boundaries( boundaries = list( mesh_points.values() ), labels = list( mesh_points.keys() ) )

for _ in range(3) :
    Th = pyff.adaptmesh( Th, iso = 1, hmax = 0.1 )

#########################
#
# Plot mesh
#
#########################

fig_mesh = figure()
ax_mesh = gca()

Th.plot_triangles( color = 'grey', alpha = .1, lw = .5 )
Th.plot_boundaries( )

legend( loc = 'center')


#########################
#
# Wave field
#
#########################

k = 2

z = Th.x + 1j*Th.y
phi = exp( -1j*k*z )

fig_phi = figure()
ax_phi = gca()
ax_phi.tricontourf( Th, real(phi) )

#########################
#
# gradient matrices
#
#########################

M = pyff.gradient_matrices( Th )

print( shape(M['grad_x']) )
print( len(Th.triangles), len(Th.x) )

proj = pyff.get_projector( Th, 'P1', 'P0' )
print(shape(proj))

dx_phi = M['grad_x']@( phi )
dy_phi = M['grad_y']@( phi )

x_P0 = proj@Th.x
y_P0 = proj@Th.y
phi_P0 = proj@phi


# quiver( x_P0, y_P0, real( dx_phi ), real( dy_phi ) )
quiver( x_P0, y_P0, imag( phi_P0.conj()*dx_phi ), imag( phi_P0.conj()*dy_phi ) )

#########################
#
# Show figures
#
#########################

for ax in [ax_mesh, ax_phi] :
    ax.axis('scaled')
    ax.axis('off')
    ax.set_xticks([])
    ax.set_yticks([])

ax_phi.set_xlim( array( ax_phi.get_xlim() )/2 )


# fig_mesh.savefig( '../../figures/gradient_wave_and_phase_mesh.svg', bbox_inches = 'tight' )
# fig_phi.savefig( '../../figures/gradient_wave_and_phase_energy_flux.svg', bbox_inches = 'tight' )



show()
