import sys
sys.path.append('./../../')

from pylab import *
import pyFreeFem as pyff

##################
#
# MESH
#
##################

fig_tau = figure()

err = 5e-3
W = 5

for mu in (.5, 1, 3):

    x = array( [0] + list( logspace( -4, log10( W/2 ), 50 ) ) )
    y = exp( -mu*x ) + exp( mu*( x - W ) ) - 1.

    x = list(x) + [W/2]
    y = list(y) + [0]

    Th = pyff.TriMesh(x,y)
    Th.add_boundary_edges( range(len(x)-1) ) # bottom
    Th.add_boundary_edges( [ len(x)-1, 0 ] ) # free surface
    Th.add_boundary_edges( [ len(x)-2, len(x)-1 ] ) # symmetry axis

    initial_Th = Th


    script = pyff.InputScript( Th = Th, err = err )
    script += '''
    fespace Vh( Th, P2 );
    fespace Vhd( Th, P1 );

    Vh u, ut;
    Vhd tauz, uOutput;

    problem Poisson( u, ut ) =
        int2d(Th)( dx(u)*dx(ut) + dy(u)*dy(ut) )
        - int2d(Th)( ut )
        + on( 1, u = 0 );

    Poisson;

    for( int i; i < 4; i++ )
    {
        Th = adaptmesh( Th, u, err = err, iso = 1 );
        Poisson;
    }

    uOutput = u;
    tauz = dy(u);

    '''
    script += pyff.OutputScript( uOutput = 'vector', tauz = 'vector', Th = 'mesh' )

    ff_output = script.get_output()
    Th = ff_output['Th']
    u = ff_output['uOutput']
    tau_z = ff_output['tauz']

    boundary_label = 1
    boundary_segments = Th.get_boundaries()[boundary_label]

    for nodes in  boundary_segments  :
        plot( Th.x[nodes][1:], tau_z[nodes][1:], label = r'$\mu=$' + str(mu) )

legend()
xscale('log'); yscale('log')
xlabel(r'Crosswise coordinate  $y$')
ylabel(r'Vertical shear stress  $\tau_z$')

####################
#
# Mesh
#
####################

for Th in ( Th, initial_Th) :
    figure()

    try :
        tricontourf( Th, u )
    except :
        Th.plot_triangles( lw = .5, alpha = .3, color = 'k' )

    Th.plot_boundaries( clip_on = False )
    legend(title = 'Boundary label', loc = 'lower left')
    axis('equal')
    xticks([])
    yticks([])
    axis('off')

fig_path = './../../figures/'

# for i in plt.get_fignums():
#     figure(i)
#     savefig(fig_path + __file__.split('/')[-1].split('.')[0] + '_' + str(i) + '.svg', bbox_inches = 'tight')

show()
