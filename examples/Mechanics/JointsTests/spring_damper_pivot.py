
from __future__ import print_function
import os,sys
import numpy as np
import math

np.set_printoptions(precision=3)

from siconos.mechanics.collision.tools import Contactor
from siconos.mechanics.joints import cast_PivotJointR
from siconos.io.mechanics_io import Hdf5
from siconos.kernel import SiconosVector, BlockVector, changeFrameAbsToBody

# An example of applying force to the axis of a joint, and applying
# spring and virtual damping by measuring position and velocity along
# the same axis.

# Note: This example is to demonstrate external measurement of joint
# positions and application of forces to dynamical systems attached to
# joints.  In practice it is better to use internal forces (fInt,
# mInt) to model joint spring-dampers, see folder
# JointsTestsWithInternalForces, and extra Relations with associated
# Non-Smooth Laws to model non-linearities such as joint stops and
# friction, see JointsTestsWithContactDetection.

# Creation of the hdf5 file for input/output
with Hdf5() as io:

    # Definition of two bars connected by a prismatic joint
    io.addPrimitiveShape('Bar', 'Box', (1, 0.1, 0.1))
    io.addObject('bar1', [Contactor('Bar')], [0,0,2],
                 orientation=[(0,0,1),np.pi/2], mass=1.0, velocity=[0,0,0,0,0,1])
    io.addObject('bar2', [Contactor('Bar',relative_translation=[0, 0.1,0]),
                          Contactor('Bar',relative_translation=[0,-0.1,0])],
                 [0,0,2],
                 orientation=[(0,0,1),np.pi/2], mass=1.0)
    io.addJoint('joint1', 'bar1', 'bar2', [0,0,0], [0,1,0], 'PivotJointR',
                absolute=False)

    # Definition of the ground
    io.addPrimitiveShape('Ground', 'Box', (5, 5, 0.1))
    io.addObject('ground', [Contactor('Ground')], [0,0,-0.05])
    io.addNewtonImpactFrictionNSL('contact', mu=0.3, e=0.0)

class Ctrl(object):
    def initialize(self, io):
        self.count = 0
        self.topo = io._model.nonSmoothDynamicalSystem().topology()
        self.ds1 = self.topo.getDynamicalSystem('bar1')
        self.ds2 = self.topo.getDynamicalSystem('bar2')
        self.joint1 = cast_PivotJointR(
            self.topo.getInteraction('joint1').relation())

        self.ds1.setIsMextExpressedInInertialFrame(True)
        self.ds2.setIsMextExpressedInInertialFrame(True)

        # Apply initial forces
        self.step()

    def step(self):
        self.count += 1

        # Make a temporary BlockVector containing both qs
        bv = BlockVector(self.ds1.q(), self.ds2.q())

        torque1 = np.zeros(3)
        torque2 = np.zeros(3)

        # Get the position and use it to project a torque vector
        # onto the DoF (spring torque)
        angle = SiconosVector(1)
        self.joint1.computehDoF(0, bv, angle, 0)

        setpoint = np.pi/4
        ang_diff = setpoint - angle.getValue(0)
        spring_torque = np.array(self.joint1.normalDoF(bv, 0)) * ang_diff * 500.0

        # Get the velocity of each body projected onto the DoF and
        # calculate their difference (damping torque)
        vel1 = self.joint1.projectVectorDoF(self.ds1.angularVelocity(True), bv, 0)
        vel2 = self.joint1.projectVectorDoF(self.ds2.angularVelocity(True), bv, 0)
        vel_diff = vel1 - vel2
        damping_torque = vel_diff * 5.0

        # Calculate total torques for each body
        torque1 += - (spring_torque + damping_torque)/2
        torque2 += + (spring_torque + damping_torque)/2

        print('applying spring-damper torques', torque1, torque2)

        self.ds1.setMExtPtr(torque1)
        self.ds2.setMExtPtr(torque2)

# Load and run the simulation
with Hdf5(mode='r+') as io:
    io.run(t0=0,
           T=5,
           h=0.001,
           theta=0.5,
           Newton_max_iter=1,
           controller=Ctrl(),
           itermax=1000,
           tolerance=1e-12,
           output_frequency=1)
