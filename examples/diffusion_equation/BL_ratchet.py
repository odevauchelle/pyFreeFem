from pylab import *
from scipy.sparse.linalg import eigs

import sys
sys.path.append('./../../')

import pyFreeFem as pyff

# create initial mesh

script = pyff.edpScript('mesh Th = square( 50, 50 );')
script += pyff.OutputScript( Th = 'mesh' )
Th = script.get_output()['Th']

Th.x -= .5
Th.y -= .5

# define FE matrices

script = pyff.InputScript( Th = 'mesh' )

script += '''
fespace Vh( Th, P2 );
Vh u, v;

real l = .2;
Vh D = 2 + tanh( x/l );
Vh beta = pow( D, -2/3. );
'''

# script += pyff.InputScript( D = 'vector', beta = 'vector' )

script += pyff.VarfScript(
    Gramian = 'int2d(Th)( u*v )',
    force = 'int2d(Th)( dy(v)*D*beta*u )',
    stiffness_classical = 'int2d(Th)( dx(u)*D*dx(v) +  dy(u)*D*dy(v) )',
    stiffness_exotic = 'int2d(Th)( dx(u)*D*dx(v) +  dy(u)*D*dy(v) + dx(D)*u*dx(v) +  dy(D)*u*dy(v) )',
    )

script += pyff.OutputScript( D = 'vector', beta = 'vector' )


# X = script.get_output( Th = Th, D = D, beta = D )['X']

# Define diffusivity


f = 0.3

p = pyff.get_projector( Th, 'P2', 'P1' )
# print(shape(p))

for diffusion_type in ('classical','exotic') :

    fig, axs = subplots( 4, 4, figsize = (12,10) )
    fig.suptitle( diffusion_type.capitalize() + ' diffusion' )
    axs = axs.flatten()
    for ax in axs :
        ax.axis('equal')
        ax.axis('off')



    # get FE matrices

    M = script.get_output( Th = Th )

    axs[0].tricontourf( Th, p@M['D'] )
    axs[0].set_title('Diffusivity', color = 'tab:red')

    axs[1].tricontourf( Th, 1./( p@M['beta'] ) )
    axs[1].set_title('Temperature', color = 'tab:red')
    
    # solve eigenvalue problem

    eigenvalues, eigenvectors = eigs( M['stiffness_' + diffusion_type] + f*M['force'], len(axs), M['Gramian'], sigma = 0 )

    print(eigenvalues)
    print(shape(eigenvectors), len(Th.x))

    for i in range( len(axs) - 2 ) :
        axs[i+2].tricontourf( Th, p@real( eigenvectors[:,i] ) )
        axs[i+2].set_title( r'$j\omega=$' + str( round( eigenvalues[i], 4 ) ) )

    # fig.savefig( diffusion_type + '_diffusion.svg', bbox_inches = 'tight' )

show()
