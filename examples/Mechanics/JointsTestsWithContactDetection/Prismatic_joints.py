#!/usr/bin/env python
from siconos.mechanics.collision.tools import Contactor
from siconos.io.mechanics_io import Hdf5
import siconos.numerics as Numerics
import siconos.kernel as Kernel

with Hdf5() as io:
    io.addPrimitiveShape('Cube', 'Box', (1, 1, 1))
    io.addPrimitiveShape('Ground', 'Box', (10, 10, .5))
    io.addNewtonImpactFrictionNSL('contact', e=0.0, mu=0.2)
    io.addObject('ground', [Contactor('Ground')], translation=[0, 0, -0.25])

    io.addObject('cube1', [Contactor('Cube')], translation=[-1, 0, 0.5], mass=0.1)
    io.addObject('cube2', [Contactor('Cube')], translation=[ 1, 0, 0.5], mass=0.1,
                 velocity=[5,0,5,0,0,0])
    # Try velocity=[5,0,0,0,0,0]) = correct behaviour
    # Try velocity=[0,0,5,0,0,0]) = correct behaviour
    # Try velocity=[0,5,0,0,0,0]) = incorrect behaviour

    io.addJoint('joint1', 'cube1', 'cube2', [0,0,1], [], 'PrismaticJointR')

with Hdf5(mode='r+') as io:
    io.run(t0=0,
           T=3,
           h=0.001,
           theta=0.50001,
           Newton_max_iter=20,
           solver=Numerics.SICONOS_FRICTION_3D_NSGS,
           itermax=10000,
           tolerance=1e-14)

    
