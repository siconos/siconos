#!/usr/bin/env python

from siconos.mechanics.collision.tools import Contactor
from siconos.io.mechanics_io import Hdf5
import siconos.numerics as Numerics
import siconos.kernel as Kernel
import numpy as np

with Hdf5() as io:
    io.addPrimitiveShape('Box', 'Box', (1, 1, 1))
    io.addPrimitiveShape('Ground', 'Box', (10, 30, 0.5))
    io.addNewtonImpactFrictionNSL('contact', e=0.5, mu=0.5)
    io.addNewtonImpactNSL('stop', e=0.8) # not a contact NSL by default (group1,2==-1)
    io.addObject('ground', [Contactor('Ground')], translation=[0, 0, -0.25])

    # We enable self-collide on the objects and their common joint to
    # show how contact behaviour and stops can interact.

    io.addObject('cube1', [Contactor('Box')], translation=[-1, 5, 2.5],
                 mass=0.1, velocity=[0,0,20,0,0,0],
                 allow_self_collide=True)

    io.addObject('cube2', [Contactor('Box')], translation=[-1, 5, 4.0],
                 mass=0.1, velocity=[0,0,0,0,0,0],
                 allow_self_collide=True)

    # We add two stops to each joint, one positive and the other
    # negative The initial velocity causes the cube to fly upwards,
    # hit the stop, and then bounce back downwards and hit the other
    # stop.  The stops are positioned at 2.0 and -2.0 and are on axis
    # 0 (the only DoF for prismatic joint.)  The prismatic vectors are
    # diagonal for joint1, and vertical for joint2, to demonstrate the
    # joint dynamics.

    io.addJoint('joint1', 'cube1', None, [0,0.707,0.707], [], 'PrismaticJointR',
                nslaws='stop', stops=[[0, 2.0, 1], [0, -2.0, -1]])

    io.addJoint('joint2', 'cube1', 'cube2', [0,0,1], [], 'PrismaticJointR',
                nslaws='stop', stops=[[0, 2.0, 1], [0, -2.0, -1]],
                allow_self_collide=True)

with Hdf5(mode='r+') as io:
    io.run(t0=0,
           T=10,
           h=0.001,
           theta=0.50001,
           Newton_max_iter=2,
           solver=Numerics.SICONOS_FRICTION_3D_NSGS,
           itermax=10000,
           tolerance=1e-14,
           # projection_itermax=3,
           # projection_tolerance=1e-5,
           # projection_tolerance_unilateral=1e-5,
           # time_stepping=Kernel.TimeSteppingDirectProjection,
           # osi=Kernel.MoreauJeanDirectProjectionOSI,
    )
