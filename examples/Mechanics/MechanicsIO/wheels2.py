#!/usr/bin/env python

#
# Example of a small driving 4-wheeled car with no steering
#

from siconos.mechanics.collision.tools import Contactor
from siconos.io.mechanics_io import Hdf5
import siconos.numerics as Numerics
import siconos.kernel as Kernel
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

    heightmap = (np.random.uniform(0,0.3,size=(100,200))
                 + [np.linspace(0,2,200)]*100).T

    # The heightmap has tangential extents 50 by 50, and its height is
    # not scaled; actual height must be reflected by values in
    # heightmap data.
    io.addHeightMap('Terrain', heightmap, (100, 50))

    # the ground object made with the ground shape. As the mass is
    # not given, it is a static object only involved in contact
    # detection.
    io.addObject('ground', [Contactor('Terrain')], translation=[40, 0, 0])

    io.addObject('body',   [Contactor('Body')], translation=[0,0,3], mass=1)

    io.addObject('wheel5', [Contactor('Wheel', relative_translation=[0,-3,0]),
                            Contactor('Wheel', relative_translation=[0,3,0])],
                 translation=[-3,0,3], mass=1)

    io.addObject('wheel6', [Contactor('Wheel', relative_translation=[0,-3,0]),
                            Contactor('Wheel', relative_translation=[0,3,0])],
                 translation=[3,0,3], mass=1)

    io.addJoint('joint6', 'wheel5', 'body', [0, 0, 0], [0, 1, 0], 'PivotJointR')
    io.addJoint('joint7', 'wheel6', 'body', [0, 0, 0], [0, 1, 0], 'PivotJointR')

class roll(object):
    def initialize(self, io):
        ang_force = 6.0
        self.io = io
        topo = io._model.nonSmoothDynamicalSystem().topology()
        self.wheels = [topo.getDynamicalSystem('wheel%d'%i) for i in [5]]
        self.wheel_const = np.array(self.wheels[0].fExt())
        self.wheel_force = Kernel.SiconosVector(3)
        self.wheel_force.setVector(0, self.wheel_const)
        self.wheel_torque = Kernel.SiconosVector(3)
        self.wheel_torque.setValue(1, ang_force)
        [w.setFExtPtr(self.wheel_force) for w in self.wheels]
        [w.setMExtPtr(self.wheel_torque) for w in self.wheels]
    def step(self):
        pass

# Run the simulation from the inputs previously defined and add
# results to the hdf5 file. The visualisation of the output may be done
# with the vview command.
with Hdf5(mode='r+') as io:

    # By default earth gravity is applied and the units are those
    # of the International System of Units.
    io.run(with_timer=False,
           controller=roll(),
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
