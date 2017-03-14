#!/usr/bin/env python

#
# Example of two cubes, one with a convex shape, one with a primitive
# shape.
#

from siconos.mechanics.collision.tools import Contactor
from siconos.io.mechanics_io import Hdf5
import siconos.numerics as Numerics

import random
import math

n_cube=3
n_row=2
n_col=2
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
    # io.addPrimitiveShape('Ground', 'Box', (200, 200, .5))
    
    io.addPrimitiveShape('MovingGround', 'Box', (200, 200, .5))
    

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
    # io.addObject('ground', [Contactor('Ground')],
    #              translation=[50, 50, 0])
    
    io.addObject('movingground', [Contactor('MovingGround')],
                 translation=[50, 50, 0],
                 mass=1.0)

    io.addBoundaryCondition('vibration', 'movingground', indices=[0,1,2,3,4,5], bc_class='HarmonicBC',
                            a=[0.0,0.0,0.0,0.0,0.0,0.0], b=[0.0,0.0,10.0,0.0,0.0,0.0],
                            omega= [0.0,0.0,math.pi,0.0,0.0,0.0],
                            phi=[0.0,0.0,0.0,0.0,0.0,0.0])
    #io.addBoundaryCondition('fixed', 'movingground', indices=[0,1,3,4,5])

    # io.addObject('left', [Contactor('Left')],
    #              translation=[0, 50., 25.])
    # io.addObject('right', [Contactor('Right')],
    #              translation=[0, -50., 25.])
    # io.addObject('rear00', [Contactor('Rear0')],
    #              translation=[25., 0., 250.])



# Run the simulation from the inputs previously defined and add
# results to the hdf5 file. The visualisation of the output may be done
# with the vview command.


nstep=100000
step=0.0005
with Hdf5(mode='r+') as io:

    # By default earth gravity is applied and the units are those
    # of the International System of Units.
    # Because of fixed collision margins used in the collision detection,
    # sizes of small objects may need to be expressed in cm or mm.
    io.run(with_timer=False,
           time_stepping=None,
           space_filter=None,
           body_class=None,
           shape_class=None,
           face_class=None,
           edge_class=None,
           gravity_scale=1,
           t0=0,
           T=nstep*step,
           h=step,
           multipoints_iterations=True,
           theta=0.50001,
           Newton_max_iter=1,
           set_external_forces=None,
           solver=Numerics.SICONOS_FRICTION_3D_NSGS,
           itermax=100,
           tolerance=1e-4,
           numerics_verbose=False,
           output_frequency=100)
