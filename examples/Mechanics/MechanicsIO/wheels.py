#!/usr/bin/env python

#
# Example of a small driving 4-wheeled car with no steering
#

from siconos.mechanics.collision.tools import Contactor
from siconos.io.mechanics_io import Hdf5
import siconos.numerics as Numerics
import numpy as np

import pydoc
# Creation of the hdf5 file for input/output
with Hdf5() as io:

    # Definition of a non smooth law. We put the objects and heightmap
    # into different collision groups so that there are object-terrain
    # collisions but no object-object collisions.
    io.addNewtonImpactFrictionNSL('contact', mu = 0.3, e = 0.0)

    # Definition of a wheel
    io.addPrimitiveShape('Wheel', 'Cylinder', (1, 0.5))

    # Definition of the car body
    io.addPrimitiveShape('Body', 'Box', (6, 5.4, 0.5))

    # Definition of the ground shape
    io.addPrimitiveShape('Ground', 'Box', (50, 50, 1))

    # the ground object made with the ground shape. As the mass is
    # not given, it is a static object only involved in contact
    # detection.
    io.addObject('ground', [Contactor('Ground')], translation=[0, 0, 0])

    io.addObject('wheel1', [Contactor('Wheel')], translation=[-3,3,2], mass=1)
    io.addObject('wheel2', [Contactor('Wheel')], translation=[3,3,2], mass=1)
    io.addObject('wheel3', [Contactor('Wheel')], translation=[-3,-3,2], mass=1)
    io.addObject('wheel4', [Contactor('Wheel')], translation=[3,-3,2], mass=1)
    io.addObject('body',   [Contactor('Body')], translation=[0,0,2], mass=1)

    io.addJoint('joint1', 'wheel1', 'body', [0, 0, 0], [0, 1, 0], 'PivotJointR')
    io.addJoint('joint2', 'wheel2', 'body', [0, 0, 0], [0, 1, 0], 'PivotJointR')
    io.addJoint('joint3', 'wheel3', 'body', [0, 0, 0], [0, 1, 0], 'PivotJointR')
    io.addJoint('joint4', 'wheel4', 'body', [0, 0, 0], [0, 1, 0], 'PivotJointR')

    # Constant value on Y-axis angular velocity (index 4) to make the wheels spin
    ang_vel = 0.0001
    io.addBoundaryCondition('spin1', 'wheel1', indices=[4],
                            bc_class='BoundaryCondition', v = [ang_vel])
    io.addBoundaryCondition('spin2', 'wheel2', indices=[4],
                            bc_class='BoundaryCondition', v = [ang_vel])
    io.addBoundaryCondition('spin3', 'wheel3', indices=[4],
                            bc_class='BoundaryCondition', v = [ang_vel])
    io.addBoundaryCondition('spin4', 'wheel4', indices=[4],
                            bc_class='BoundaryCondition', v = [ang_vel])

# Run the simulation from the inputs previously defined and add
# results to the hdf5 file. The visualisation of the output may be done
# with the vview command.
with Hdf5(mode='r+') as io:

    # By default earth gravity is applied and the units are those
    # of the International System of Units.
    io.run(with_timer=False,
           t0=0,
           T=20,
           h=0.001,
           theta=0.50001,
           Newton_max_iter=1,
           solver=Numerics.SICONOS_FRICTION_3D_NSGS,
           itermax=10000,
           tolerance=1e-8,
           numerics_verbose=False,
           output_frequency=None)
