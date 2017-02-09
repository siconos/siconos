#!/usr/bin/env python

#
# Example of a small driving 4-wheeled car with no steering
#

from siconos.mechanics.collision.tools import Contactor
from siconos.io.mechanics_io import Hdf5
import siconos.numerics as Numerics
import siconos.kernel as Kernel
import numpy as np

bumpy_terrain = True  # false = flat ground
use_torque = True     # false = boundary condition

# Creation of the hdf5 file for input/output
with Hdf5() as io:

    # Definition of a non smooth law. We put the objects and heightmap
    # into different collision groups so that there are object-terrain
    # collisions but no object-object collisions.
    io.addNewtonImpactFrictionNSL('contact', mu = 0.6, e = 0.0)

    # Definition of a wheel (Radius = 1.5, width = 0.5)
    io.addPrimitiveShape('Wheel', 'Cylinder', (1.5, 0.5))

    # Definition of the car body
    io.addPrimitiveShape('Body', 'Box', (6, 5.4, 0.5))

    if not bumpy_terrain:
        # Definition of the ground shape
        io.addPrimitiveShape('Ground', 'Box', (50, 50, 1))

        # The ground object made with the ground shape. As the mass is
        # not given, it is a static object only involved in contact
        # detection.
        io.addObject('ground', [Contactor('Ground')], translation=[0, 0, 0])

    else:
        # The heightmap is a 100x200 grid with a height of 1.
        heightmap = (np.random.uniform(0,0.3,size=(100,200))
                     + [np.linspace(0,2,200)]*100
                     + [np.sin(np.linspace(0,np.pi*12,200))]*100).T

        # It has tangential extents 100 by 50.
        io.addHeightMap('Terrain', heightmap, (100, 50))

        # It is also a static object (no mass)
        io.addObject('ground', [Contactor('Terrain')], translation=[40, 0, 0])

    h = 3
    io.addObject('wheel1', [Contactor('Wheel')], translation=[-3, 3, h], mass=1)
    io.addObject('wheel2', [Contactor('Wheel')], translation=[ 3, 3, h], mass=1)
    io.addObject('wheel3', [Contactor('Wheel')], translation=[-3,-3, h], mass=1)
    io.addObject('wheel4', [Contactor('Wheel')], translation=[ 3,-3, h], mass=1)
    io.addObject('body',   [Contactor('Body')],  translation=[ 0, 0, h], mass=1)

    io.addJoint('joint1', 'wheel1', 'body', [0, 0, 0], [0, 1, 0], 'PivotJointR')
    io.addJoint('joint2', 'wheel2', 'body', [0, 0, 0], [0, 1, 0], 'PivotJointR')
    io.addJoint('joint3', 'wheel3', 'body', [0, 0, 0], [0, 1, 0], 'PivotJointR')
    io.addJoint('joint4', 'wheel4', 'body', [0, 0, 0], [0, 1, 0], 'PivotJointR')

    # Constant value on Y-axis angular velocity (index 4) to make the wheels spin
    if not use_torque:
        ang_vel = 8.0
        io.addBoundaryCondition('spin1', 'wheel1', indices=[4],
                                bc_class='BoundaryCondition', v = [ang_vel])
        io.addBoundaryCondition('spin2', 'wheel2', indices=[4],
                                bc_class='BoundaryCondition', v = [ang_vel])
        io.addBoundaryCondition('spin3', 'wheel3', indices=[4],
                                bc_class='BoundaryCondition', v = [ang_vel])
        io.addBoundaryCondition('spin4', 'wheel4', indices=[4],
                                bc_class='BoundaryCondition', v = [ang_vel])

# A controller to give constant torque to the wheels
class roll(object):
    def initialize(self, io):
        ang_force = 6.0
        self.io = io
        topo = io._model.nonSmoothDynamicalSystem().topology()
        self.wheels = [topo.getDynamicalSystem('wheel%d'%i) for i in [1,2,3,4]]
        self.wheel_const = np.array(self.wheels[0].fExt())
        self.wheel_force = Kernel.SiconosVector(3)
        self.wheel_force.setVector(0, self.wheel_const)
        self.wheel_torque = Kernel.SiconosVector(3)
        self.wheel_torque.setValue(1, ang_force)

        # Same force and torque to all wheels
        [w.setFExtPtr(self.wheel_force) for w in self.wheels]
        [w.setMExtPtr(self.wheel_torque) for w in self.wheels]
    def step(self):
        pass

controller = None
if use_torque:
    controller = roll()

# Run the simulation from the inputs previously defined and add
# results to the hdf5 file. The visualisation of the output may be done
# with the vview command.
with Hdf5(mode='r+') as io:

    # By default earth gravity is applied and the units are those
    # of the International System of Units.
    io.run(with_timer=False,
           t0=0,
           T=20,
           h=0.005,
           controller=controller,
           theta=0.50001,
           Newton_max_iter=1,
           solver=Numerics.SICONOS_FRICTION_3D_NSGS,
           itermax=1000,
           tolerance=1e-4,
           numerics_verbose=False,
           output_frequency=None)
