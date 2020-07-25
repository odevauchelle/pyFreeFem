from pylab import *

import sys
sys.path.append('./../../')

import pyFreeFem as pyff

# Mesh

nbPts = 25

script = pyff.InputScript( nbPts = nbPts )
script += pyff.edpScript('mesh Th = square( nbPts, nbPts );')
script += pyff.OutputScript( Th = 'mesh' )

# Finite element space and matrices

script += pyff.edpScript('fespace Vh( Th, P1 );')
script += pyff.VarfScript(
    stiffness = 'int2d(Th)( dx(u)*dx(v) +  dy(u)*dy(v) )',
    Gramian = 'int2d(Th)( u*v )',
    BoundaryGramian = 'int1d(Th, 1, 2, 3, 4 )( u*v )'
    )

ff_out = script.get_output()

# Eigenvalue problem
from scipy.sparse.linalg import eigs

epsilon = 1e-6

eigenvalues, eigenvectors = eigs( ff_out['stiffness'] - 1/epsilon*ff_out['BoundaryGramian'], 9, ff_out['Gramian'], sigma = 5 )

# Plot

ncols = 3

Th = ff_out['Th']

figure( figsize = (6,6) )

Th.plot_triangles(color = 'k', lw = .5, alpha = .2 )
Th.plot_boundaries()
legend( title = 'Boundaries', ncol = 1, loc='center', bbox_to_anchor=(0.2, 0.2), framealpha = 1 )
axis('equal')
xticks([]); yticks([])
axis('off')

savefig( '../../figures/' + __file__.split('/')[-1].split('.')[0] + '_mesh.svg' , bbox_inches = 'tight' )



fig, axs = subplots( nrows = int( len( eigenvalues )/ncols ), ncols = ncols, figsize = (7,7) )
axs = axs.flatten()


for i, eigenvalue in enumerate( eigenvalues ) :

    omega = sqrt( real( eigenvalue ) )

    ax = axs[i]
    ax.tricontourf( Th, real( eigenvectors[:,i] ) )
    ax.set_title(r'$\omega=$' + str( round( omega, 2 ) ) )

for ax in axs :

    Th.plot_boundaries( color = 'grey', ax = ax )
    ax.axis('equal')
    ax.axis('off')
    ax.set_xticks([])
    ax.set_yticks([])


savefig( '../../figures/' + __file__.split('/')[-1].split('.')[0] + '.svg' , bbox_inches = 'tight' )
show()
