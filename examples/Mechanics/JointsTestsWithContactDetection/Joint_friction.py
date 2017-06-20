#!/usr/bin/env python

from siconos.mechanics.collision.tools import Contactor
from siconos.io.mechanics_io import Hdf5
import siconos.numerics as Numerics
import siconos.kernel as Kernel
from siconos.mechanics.collision.bullet import SiconosBulletOptions
import numpy as np

with Hdf5() as io:
    io.addPrimitiveShape('BigBox', 'Box', (1, 1, 1))
    io.addPrimitiveShape('SmallBox', 'Box', (0.3,0.3,0.3))
    io.addNewtonImpactFrictionNSL('contact', e=0.0, mu=0.0)

    io.addObject('heavy1', [Contactor('SmallBox')], translation=[0, 0, 4], mass=1.0,
                 velocity=[0,0,15,0,0,0])
    io.addObject('heavy2', [Contactor('SmallBox')], translation=[0, 0, -4], mass=1.0,
                 velocity=[0,0,15,0,0,0])

    io.addObject('big', [Contactor('BigBox')], translation=[0, 0, 0], mass=0.1,
                 velocity=[0,0,0,0,0,0])

    # We create a prismatic joint with friction on its axis.  In fact
    # as there is no "normal force" in a 1-D friction, it is expressed
    # as a lower- and upper-bound of permissible velocity by means of
    # a RelayNSL, here between [-0.1, 0.1].

    io.addRelayNSL('friction', lb=-0.1, ub=0.1)

    io.addJoint('joint1', 'big', None, [0,0,1], [], 'PrismaticJointR',
                friction='friction')

with Hdf5(mode='r+') as io:
    io.run(t0=0,
           T=8,
           h=0.001,
           theta=0.50001,
           Newton_max_iter=2,
           solver=Numerics.SICONOS_FRICTION_3D_NSGS,
           itermax=1000,
           tolerance=1e-4,
    )
