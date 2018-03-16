#!/usr/bin/env python

#
# Example of two cubes, one with a convex shape, one with a primitive
# shape.
#

from siconos.mechanics.collision.tools import Contactor
from siconos.io.mechanics_run import MechanicsHdf5Runner
import siconos.numerics as Numerics

import random
import math

n_cube=3
n_row=2
n_col=2
# Creation of the hdf5 file for input/output
with MechanicsHdf5Runner() as io:
    for i in range(n_row):
        for j in range(n_col):
            for n in range(n_cube):
            # Definition of a cube as a convex shape
                io.add_convex_shape('CubeCS'+str(n)+'_'+str(i)+'_'+str(j), [ (-1.0, 1.0, -1.0),
                                                                           (-1.0, -1.0, -1.0),
                                                                           (-1.0, -1.0, 1.0),
                                                                           (-1.0, 1.0, 1.0),
                                                                           (1.0, 1.0, 1.0),
                                                                           (1.0, 1.0, -1.0),
                                                                           (1.0, -1.0, -1.0),
                                                                           (1.0, -1.0, 1.0)])



    # Alternative to the previous convex shape definition.
    #io.add_primitive_shape('CubePrim', 'Box', (2, 2, 2))

    # Definition of the ground shape
    # io.add_primitive_shape('Ground', 'Box', (200, 200, .5))
    
    io.add_primitive_shape('MovingGround', 'Box', (200, 200, .5))
    

    # Definition of the left shape
    # io.add_primitive_shape('Left', 'Box', (100, 0.5, 50.))

    # Definition of the right shape
    #io.add_primitive_shape('Right', 'Box', (100, 0.5, 50.))

    # Definition of the rear shape
    #io.add_primitive_shape('Rear0', 'Box', (0.5, 100., 50.))

    # Definition of the front shape
    #io.add_primitive_shape('Front', 'Box', (100, 0.5, 50.))

    # Definition of a non smooth law. As no group ids are specified it
    # is between contactors of group id 0.
    io.add_Newton_impact_friction_nsl('contact', mu=0.3)

    # The cube object made with an unique Contactor : the cube shape.
    # As a mass is given, it is a dynamic system involved in contact
    # detection and in the simulation.  With no group id specified the
    # Contactor belongs to group 0
    for i in range(n_row):
        for j in range(n_col):
            for n in range(n_cube):
                io.add_object('cubeCS'+str(n)+'_'+str(i)+'_'+str(j), [Contactor('CubeCS'+str(n)+'_'+str(i)+'_'+str(j))],
                             translation=[3.0*i, 3.0*j, 2.05*(n+1)],
                             velocity=[10*(1.0+2.0*(random.random()-1.0)/2.0), 10*(1.0+2.0*(random.random()-1.0)/2.0), 0, 1, 1, 1],
                             mass=1)

    # io.add_object('cube2', [Contactor('CubePrim')], translation=[0, 3, 2],
    #              velocity=[10, 0, 0, 1, 1, 1],
    #              mass=1)

    # the ground object made with the ground shape. As the mass is
    # not given, it is a static object only involved in contact
    # detection.
    # io.add_object('ground', [Contactor('Ground')],
    #              translation=[50, 50, 0])
    
    io.add_object('movingground', [Contactor('MovingGround')],
                 translation=[50, 50, 0],
                 mass=1.0)

    io.add_boundary_condition('vibration', 'movingground', indices=[0,1,2,3,4,5], bc_class='HarmonicBC',
                            a=[0.0,0.0,0.0,0.0,0.0,0.0], b=[0.0,0.0,10.0,0.0,0.0,0.0],
                            omega= [0.0,0.0,math.pi,0.0,0.0,0.0],
                            phi=[0.0,0.0,0.0,0.0,0.0,0.0])
    #io.add_boundary_condition('fixed', 'movingground', indices=[0,1,3,4,5])

    # io.add_object('left', [Contactor('Left')],
    #              translation=[0, 50., 25.])
    # io.add_object('right', [Contactor('Right')],
    #              translation=[0, -50., 25.])
    # io.add_object('rear00', [Contactor('Rear0')],
    #              translation=[25., 0., 250.])



# Run the simulation from the inputs previously defined and add
# results to the hdf5 file. The visualisation of the output may be done
# with the vview command.


nstep=100000
step=0.0005
with MechanicsHdf5Runner(mode='r+') as io:

    # By default earth gravity is applied and the units are those
    # of the International System of Units.
    # Because of fixed collision margins used in the collision detection,
    # sizes of small objects may need to be expressed in cm or mm.
    io.run(with_timer=False,
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
