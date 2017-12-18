#!/usr/bin/env python

#
# Example of two cubes, one with a convex shape, one with a primitive
# shape.
#

from siconos.mechanics.collision.tools import Contactor
from siconos.io.mechanics_io import Hdf5
import siconos.numerics as Numerics
import siconos.kernel as Kernel

import random

import siconos
options = siconos.mechanics.collision.bullet.SiconosBulletOptions()
options.worldScale = 1.0
options.perturbationIterations = 7
options.minimumPointsPerturbationThreshold = 7

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
    io.addPrimitiveShape('Ground', 'Box', (200, 200, .5))

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



# Run the simulation from the inputs previously defined and add
# results to the hdf5 file. The visualisation of the output may be done
# with the vview command.


nstep=2000
step=0.005
with Hdf5(mode='r+') as io:

    # By default earth gravity is applied and the units are those
    # of the International System of Units.
    # Because of fixed collision margins used in the collision detection,
    # sizes of small objects may need to be expressed in cm or mm.
    io.run(with_timer=False,
           time_stepping=Kernel.TimeSteppingDirectProjection,
           osi=Kernel.MoreauJeanDirectProjectionOSI,
           space_filter=None,
           body_class=None,
           shape_class=None,
           face_class=None,
           edge_class=None,
           gravity_scale=1,
           options = options,
           t0=0,
           T=nstep*step,
           h=step,
           theta=0.50001,
           Newton_max_iter=1,
           Newton_update_interactions=False,
           set_external_forces=None,
           solver=Numerics.SICONOS_FRICTION_3D_NSGS,
           itermax=100,
           tolerance=1e-4,
           numerics_verbose=False,
           output_frequency=1,
           projection_itermax=5,
           projection_tolerance=1e-8,
           projection_tolerance_unilateral=1e-8,
    )
