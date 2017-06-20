#!/usr/bin/env python

#
# Example of two cubes, one with a convex shape, one with a primitive
# shape.
#

from siconos.mechanics.collision.tools import Contactor
from siconos.io.mechanics_io import Hdf5
import siconos.numerics as Numerics

import random

n_cube=4
n_row=4
n_col=4
# Creation of the hdf5 file for input/output
with Hdf5() as io:
    for i in range(n_row):
        for j in range(n_col):
            for n in range(n_cube):
            # Definition of a cube as a convex shape
                io.addConvexShape('CubeCS'+str(n)+'_'+str(i)+'_'+str(j), [ (-1.0, 1.0, -1.0),
                                                                           (-1.0, -1.0, -1.0),
                                                                           (-1.0, -1.0, 1.0),
                                                                           (-1.0, 1.0, 1.0),
                                                                           (1.0, 1.0, 1.0),
                                                                           (1.0, 1.0, -1.0),
                                                                           (1.0, -1.0, -1.0),
                                                                           (1.0, -1.0, 1.0)])



    # Alternative to the previous convex shape definition.
    #io.addPrimitiveShape('CubePrim', 'Box', (2, 2, 2))

    # Definition of the ground shape
    io.addPrimitiveShape('Ground', 'Box', (200, 200, .5), insideMargin=0.04)

    # Definition of the left shape
    # io.addPrimitiveShape('Left', 'Box', (100, 0.5, 50.))

    # Definition of the right shape
    #io.addPrimitiveShape('Right', 'Box', (100, 0.5, 50.))

    # Definition of the rear shape
    #io.addPrimitiveShape('Rear0', 'Box', (0.5, 100., 50.))

    # Definition of the front shape
    #io.addPrimitiveShape('Front', 'Box', (100, 0.5, 50.))

    # Definition of a non smooth law. As no group ids are specified it
    # is between contactors of group id 0.
    io.addNewtonImpactFrictionNSL('contact', mu=0.3)

    # The cube object made with an unique Contactor : the cube shape.
    # As a mass is given, it is a dynamic system involved in contact
    # detection and in the simulation.  With no group id specified the
    # Contactor belongs to group 0
    for i in range(n_row):
        for j in range(n_col):
            for n in range(n_cube):
                io.addObject('cubeCS'+str(n)+'_'+str(i)+'_'+str(j), [Contactor('CubeCS'+str(n)+'_'+str(i)+'_'+str(j))],
                             translation=[3.0*i, 3.0*j, 2.05*(n+1)],
                             velocity=[10*(1.0+2.0*(random.random()-1.0)/2.0), 10*(1.0+2.0*(random.random()-1.0)/2.0), 0, 1, 1, 1],
                             mass=1)

    # io.addObject('cube2', [Contactor('CubePrim')], translation=[0, 3, 2],
    #              velocity=[10, 0, 0, 1, 1, 1],
    #              mass=1)

    # the ground object made with the ground shape. As the mass is
    # not given, it is a static object only involved in contact
    # detection.
    io.addObject('ground', [Contactor('Ground')],
                 translation=[50, 50, 0])
    # io.addObject('left', [Contactor('Left')],
    #              translation=[0, 50., 25.])
    # io.addObject('right', [Contactor('Right')],
    #              translation=[0, -50., 25.])
    # io.addObject('rear00', [Contactor('Rear0')],
    #              translation=[25., 0., 250.])
