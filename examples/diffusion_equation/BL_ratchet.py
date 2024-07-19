from pylab import *
from scipy.sparse.linalg import eigs

import sys
sys.path.append('./../../')

import pyFreeFem as pyff

# create initial mesh

script = pyff.edpScript('mesh Th = square( 20, 20 );')
script += pyff.OutputScript( Th = 'mesh' )
Th = script.get_output()['Th']

Th.x -= .5
Th.y -= .5

# define FE matrices

script = pyff.InputScript( Th = 'mesh' )

script += '''
fespace Vh( Th, P1 );
Vh u, v;

'''

script += pyff.InputScript( D = 'vector', beta = 'vector' )


script += pyff.VarfScript(
    Gramian = 'int2d(Th)( u*v )',
    force = 'int2d(Th)( dy(v)*D*beta*u )',
    stiffness_classical = 'int2d(Th)( dx(u)*D*dx(v) +  dy(u)*D*dy(v) )',
    stiffness_exotic = 'int2d(Th)( dx(u)*D*dx(v) +  dy(u)*D*dy(v) + dx(D)*u*dx(v) +  dy(D)*u*dy(v) )',
    )

# script += pyff.OutputScript( X = 'vector' )


# X = script.get_output( Th = Th, D = D, beta = D )['X']

# Define diffusivity

l = .2
D = 2 + tanh( Th.x/l )
beta = 1/D**(2/3)
f = 0.3

# p = pyff.get_projector( Th, 'P2', 'P1' )
# print(shape(p))

for diffusion_type in ('classical','exotic') :

    fig, axs = subplots( 2, 3 )
    fig.suptitle( diffusion_type.capitalize() + ' diffusion' )
    axs = axs.flatten()
    for ax in axs :
        ax.axis('equal')
        ax.axis('off')


    axs[0].tricontourf( Th, D )
    axs[0].set_title('Diffusivity')

    # get FE matrices

    M = script.get_output( Th = Th, D = D, beta = D )

    # solve eigenvalue problem

    eigenvalues, eigenvectors = eigs( M['stiffness_' + diffusion_type] + f*M['force'], 9, M['Gramian'], sigma = 0 )

    print(eigenvalues)
    print(shape(eigenvectors), len(Th.x))

    for i in range( len(axs) - 1 ) :
        axs[i+1].tricontourf( Th, real( eigenvectors[:,i] ) )
        axs[i+1].set_title( r'$j\omega=$' + str( round( eigenvalues[i], 4 ) ) )

    # fig.savefig( diffusion_type + '_diffusion.svg', bbox_inches = 'tight' )

show()
