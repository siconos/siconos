
import os,sys
import numpy
import math

from siconos.mechanics.collision.tools import Contactor
from siconos.io.mechanics_io import Hdf5
import siconos.io.mechanics_io
import siconos.numerics as Numerics

# A collection of box stacks for stress-testing Siconos solver with
# chains of contacts.

# Creation of the hdf5 file for input/output
with Hdf5() as io:

    width,depth,height = 1,1,1
    io.addPrimitiveShape('Box', 'Box', [width,depth,height])

    k = 0
    sep = 0.01
    def make_stack(X, Y, N, M, W):
        global k
        z = height/2.0
        while W > 0:
            for i in range(N):
                for j in range(M):
                    x = (i-N/2.0)*(width+sep) + X
                    y = (j-M/2.0)*(depth+sep) + Y
                    io.addObject('box%03d'%k, [Contactor('Box')],
                                 translation=[x,y,z],
                                 mass=1.0)
                    k += 1
            N = N - 1 if N > 1 else N
            M = M - 1 if M > 1 else M
            W = W - 1
            z += height + sep

    # A column
    make_stack(0, -10, 1, 1, 5)

    # A pyramid
    make_stack(0, 0, 5, 5, 5)

    # A wall
    make_stack(0, 10, 1, 5, 5)

    # Definition of the ground
    io.addPrimitiveShape('Ground', 'Box', (50, 50, 0.1))
    io.addObject('ground', [Contactor('Ground')], [0,0,-0.05])

    # Enable to smash the wall
    # io.addPrimitiveShape('Ball', 'Sphere', [1,])
    # io.addObject('WreckingBall', [Contactor('Ball')],
    #              translation=[30,0,3], velocity=[-30,0,2,0,0,0],
    #              mass=10)

    # Definition of a non smooth law. As no group ids are specified it
    # is between contactors of group id 0.
    io.addNewtonImpactFrictionNSL('contact', mu=0.3)

# Load and run the simulation
with Hdf5(mode='r+') as io:
    io.run(t0=0,
           T=10,
           h=0.01,
           theta=0.5,
           Newton_max_iter=1,
           solver=Numerics.SICONOS_FRICTION_3D_NSGS,
           itermax=1000,
           tolerance=1e-12,
           output_frequency=1)
