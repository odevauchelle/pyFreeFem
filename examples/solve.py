import sys
sys.path.append('./../')

import pyFreeFem as pyff
import matplotlib.pyplot as pp

script = pyff.edpScript('''
real smallRadius = .3;
border outerCircle( t = 0, 2*pi ){ x = cos(t); y = 0.8*sin(t); }
border innerCircle( t = 2*pi, 0 ){ x = .5 + smallRadius*cos(t); y = smallRadius*sin(t); }
mesh Th = buildmesh( outerCircle(100) + innerCircle(40) );

fespace Vh( Th, P1 );
Vh u,v;
''')

script += pyff.edpOutput( type = 'mesh', name = 'Th' )

matrices = {
    'stiffness' : 'int2d(Th)( dx(u)*dx(v) +  dy(u)*dy(v) )',
    'Grammian' : 'int2d(Th)( u*v )',
    'boundary_Grammian' : 'int1d(Th, 1, 2)( u*v )'
}

for matrix_name in matrices.keys() :
    script += pyff.VarfBlock( name = matrix_name, varf = matrices[matrix_name] )

ff_output = script.get_output()
Th = ff_output['Th']

Th.plot_triangles( color = 'k', alpha = .2, lw = .5 )
Th.plot_boundaries()

##########################################################################################################

pp.axis('equal')
pp.axis('off')
pp.xticks([])
pp.yticks([])

pp.savefig( '../figures/' + __file__.split('/')[-1].split('.')[0] + '.svg' , bbox_inches = 'tight' )

x_lim, y_lim = pp.gca().get_xlim(), pp.gca().get_ylim()

pp.show()

##########################################################################################################

from scipy.sparse.linalg import spsolve
import numpy as np

epsilon = 1e-4
M = - ff_output['stiffness'] + 1./epsilon*ff_output['boundary_Grammian']
Source = ff_output['Grammian']*np.array( [1]*len( Th.x ) )
u = spsolve( M, Source )

pp.tricontourf( Th, u )
Th.plot_boundaries( color = 'black' )

##########################################################################################################

pp.xlim(x_lim); pp.ylim(y_lim)
pp.axis('off')
pp.xticks([])
pp.yticks([])

pp.savefig( '../figures/' + __file__.split('/')[-1].split('.')[0] + '_2.svg' , bbox_inches = 'tight' )

pp.show()

##########################################################################################################
