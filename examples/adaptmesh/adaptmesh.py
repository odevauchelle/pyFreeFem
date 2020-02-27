from pylab import *

import sys
sys.path.append('./../../')

import pyFreeFem as pyff

# create initial mesh

script = pyff.edpScript('mesh Th = square(5, 5);')
script += pyff.edpOutput( data_type = 'mesh', name = 'Th' )
Th = script.get_output()['Th']

Th_initial = Th

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

# plots

fig_index = 0

for Th in [ Th_initial, Th ] :

    figure( fig_index, figsize = (5,5))

    Th.plot_triangles(color = 'k', lw = .5, alpha = .2)
    Th.plot_boundaries()

    axis('equal')
    axis('off')
    xticks([]), yticks([])

    # savefig( './../../figures/' + __file__.split('/')[-1].split('.')[0] + '_' + str(fig_index) + '.svg' , bbox_inches = 'tight' )

    fig_index += 1

show()
