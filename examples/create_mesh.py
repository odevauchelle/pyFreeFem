import sys
sys.path.append('./../')

import pyFreeFem as pyff
import matplotlib.pyplot as pp

edp_str = '''
border Circle( t = 0, 2*pi ){ x = cos(t); y = sin(t); }
mesh Th = buildmesh( Circle(10) );
'''

edp_str += pyff.export_mesh_edp() # adds a few lines to edp string

FreeFem_output = pyff.run_FreeFem( edp_str )

mesh = pyff.FreeFem_str_to_mesh( FreeFem_output, simple_boundaries = True )

triangles_plot = mesh.plot_triangles( labels = 'index' )
mesh.plot_nodes( labels = 'index', color = triangles_plot[0].get_color() )
mesh.plot_edges( labels = 'label', color = 'red')

for key in mesh.boundary_edges.keys() :
    print(key,mesh.boundary_edges[key] )


pp.axis('equal')
pp.axis('off')
pp.xticks([])
pp.yticks([])

pp.savefig( '../figures/' + __file__.split('/')[-1].split('.')[0] + '.svg' , bbox_inches = 'tight' )

pp.show()
