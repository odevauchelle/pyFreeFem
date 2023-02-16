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

L_box = 7.
H_box = 1.5
L_boat = 1.
boat_draft = .3

x_boat = linspace(-.5,.5,50)[::-1]*L_boat
y_boat = -boat_draft*( 1 - (2*x_boat)**2 )**2*( 1 - .0*x_boat )

mesh_points = dict( boat = array([x_boat, y_boat]).T )
mesh_points['left_surface'] = [ mesh_points['boat'][-1], [-L_box/2,0] ]
mesh_points['left_wall'] = [ mesh_points['left_surface'][-1], [-L_box/2,-H_box] ]
mesh_points['bottom'] = [ mesh_points['left_wall'][-1], [L_box/2,-H_box] ]
mesh_points['right_wall'] = [ mesh_points['bottom'][-1], [L_box/2,0] ]
mesh_points['right_surface'] = [ mesh_points['right_wall'][-1], [L_boat/2,0] ]

Th = pyff.TriMesh_from_boundaries( boundaries = list( mesh_points.values() ), labels = list( mesh_points.keys() ) )

for _ in range(3) :
    Th = pyff.adaptmesh( Th, iso = 1, hmax = 0.1 )

#########################
#
# Prepare figure
#
#########################

fig_phys, axs_phys = subplots( nrows = 3, sharex = True )

for ax in axs_phys :
    Th.plot_triangles( ax = ax, color = 'w', alpha = .1, lw = .5 )
    Th.plot_boundaries( ax = ax, labels = ['boat'], color = 'red' )

    ax.axis('scaled')
    ax.axis('off')

#########################
#
# Finite elements matrices
#
#########################


script = pyff.InputScript( Th = Th )

boundary_name2int, boundary_int2name = Th.get_boundary_label_conversion()

free_surface_int_labels = [ str( boundary_name2int[ name ] ) for name in ( 'left_surface', 'right_surface' ) ]
walls_int_labels = [ str( boundary_name2int[ name ] ) for name in ( 'left_wall', 'right_wall' ) ]

script += pyff.edpScript('fespace Vh( Th, P1 );')
script += pyff.VarfScript(
    stiffness = 'int2d(Th)( dx(u)*dx(v) +  dy(u)*dy(v) )',
    Gramian = 'int2d(Th)( u*v )',
    FreeSurfaceGramian = 'int1d(Th,' + ','.join( free_surface_int_labels ) +  ' )( u*v )',
    WallsGramian = 'int1d(Th,' + ','.join( walls_int_labels ) +  ' )( u*v )',
    BoatGramian = 'int1d(Th,' + str( boundary_name2int['boat'] ) + ' )( u*v*N.y )',
    BoatForce = 'int1d(Th,' + str( boundary_name2int['boat'] ) + ' )( u*dx(v)*N.y )',
    )

ff_out = script.get_output()

#########################
#
# solve
#
#########################

k = 5
omega = sqrt( k*tanh( k*H_box ) ) #dispersion relation

M = ff_out['stiffness'] - omega**2*ff_out['FreeSurfaceGramian'] + 1j*k*ff_out['WallsGramian']

h1_hull = dict( heave = 0*Th.x + 1, pitch = 1.*Th.x )
phi = {}

for i, name in enumerate( h1_hull.keys() ) :

    phi[name] = spsolve( M, 1j*omega*ff_out['BoatGramian']*h1_hull[name] )

    axs_phys[i].tricontourf( Th, real( phi[name] ) )
    axs_phys[i].set_title(name.capitalize())

##########################
#
# Hysteretic thrust
#
##########################

hyst_thrust = []
motions = list( phi.keys())

for motion in [ motions, motions[::-1] ] :
    hyst_thrust += [ ( 1j*omega*phi[motion[0]] ).conj().T@ff_out['BoatForce']@h1_hull[motion[1]] ]

theta_max = angle( conj( hyst_thrust[0] ) -  hyst_thrust[1] )  + pi/2
print('Optimal phase:', theta_max )

phi_max = phi['heave'] + exp(1j*theta_max)*phi['pitch']

axs_phys[-1].tricontourf( Th, real( phi_max ) )
axs_phys[-1].set_title('Best heave-pitch phase: ' + str( round( theta_max/pi, 2 ) ) + r'$\pi$' )

print('Force from hysteretic thrust:', imag( exp(1j*theta_max)*hyst_thrust[0] + exp(-1j*theta_max)*hyst_thrust[1] )/2  )

##########################
#
# Amplitude of the waves
#
##########################

phi_b = dict( left = {}, right = {} )


figure()
ax_w = gca()
ax_w.set_xlabel(r'$|\phi|$')
ax_w.set_ylabel(r'$y$')

for name in phi_b.keys() :

    phi_b[name]['nodes'] = Th.get_boundaries()[ name + '_wall']

    ax_w.plot( abs( phi_max[phi_b[name]['nodes']] ), Th.y[ phi_b[name]['nodes'] ], label = name )

    phi_b[name]['A'] = abs( k*trapz( phi_max[phi_b[name]['nodes']], Th.y[ phi_b[name]['nodes'] ] ) )

    print( name, phi_b[name]['A'] )

print('Force from wave amplitude:', k*( phi_b['left']['A']**2 - phi_b['right']['A']**2 )/2 )

ax_w.legend()
# ax_w.set_xscale('log')



##########################
#
# show figures
#
##########################
# fig_phys.savefig('gunwale_bobbing.pdf', bbox_inches = 'tight')
show()
