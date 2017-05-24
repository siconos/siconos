#!/usr/bin/env python

from __future__ import print_function

from siconos.mechanics.collision.tools import Contactor
from siconos.io.mechanics_io import Hdf5
import siconos.numerics as Numerics
import siconos.kernel as Kernel
import siconos.mechanics as Mechanics
import numpy as np

# Configuration of two stops, at 3/4 pi on either side of the initial position.
# The values are [axis, position, direction], where direction must be -1 or 1.
stops = [[0, -np.pi*3/4, 1], [0, np.pi*3/4, -1]]

# Initial rotational velocity of the body attached to the joint.
twist = 10

with Hdf5() as io:
    # self-collide property: if set to true, collisions between bodies
    # connected by joints are allowed.  In this example, the joint
    # keeps the bodies in an overlapping state, which leads to
    # ill-posed problems for the solver.
    self_collide = False

    # A "bar" connected to a "post".  A "knob" is attached to the post
    # for visual reference -- otherwise no rotation can be seen, since
    # it is cylindrical.
    io.addPrimitiveShape('Bar', 'Box', (1, 0.1, 0.1))
    io.addPrimitiveShape('Post', 'Cylinder', (0.05, 1))
    io.addPrimitiveShape('Knob', 'Box', (0.2, 0.05, 0.05))
    io.addPrimitiveShape('Ground', 'Box', (4, 4, .5))

    # Ground is defined for visual reference
    io.addObject('ground', [Contactor('Ground')],
                 translation=[0, 0, 0])

    # We define a contact law even though this simulation should not
    # feature contact.
    io.addNewtonImpactFrictionNSL('contact', e=0.7, mu=0.02)

    # This law is used to specify "bouncy" stops on joint1.
    io.addNewtonImpactNSL('stop', e=0.8)

    # Very high friction on the second joint causes "almost-fixed"
    # behaviour, resulting in a bounce when the first joint hits the
    # stop, with energy partially absorbed by a small movement of the
    # second joint.
    io.addRelayNSL('friction', lb=-3.0, ub=3.0)

    # The objects, with self-collision disabled as noted above.
    io.addObject('bar', [Contactor('Bar')], translation=[0.45, 0.45, 3], mass=10,
                 allow_self_collide = self_collide, velocity=[0,0,0,0,twist,0])
    io.addObject('post', [Contactor('Post'),
                          Contactor('Knob', relative_translation=[0.1,0,0])],
                 translation=[0, 0, 3], mass=1,
                 allow_self_collide = self_collide)

    # Connect the two bodies by a pivot joint
    io.addJoint('joint1', 'bar', 'post', [-0.45,0,0], [0,1,0], 'PivotJointR',
                allow_self_collide = self_collide,
                nslaws='stop', stops=stops)

    # Joint from "bar" to the world reference frame, to keep things from falling.
    io.addJoint('joint2', 'post', None, [0,0,0], [0,1,0],
                'PivotJointR', absolute=False,
                friction = 'friction')

    # For fully fixed behaviour replace with a FixedJointR.
    # io.addJoint('joint2', 'post', None, [0,0,0], [], 'FixedJointR')

# We define a "controller" here to show how to measure the angle of
# the joint using computehDoF.
class Ctrl(object):
    def initialize(self, io):
        self.nsds = io._model.nonSmoothDynamicalSystem()
        self.topo = self.nsds.topology()
        self.joint1_inter = self.topo.getInteraction('joint1')
        self.joint1 = Mechanics.joints.cast_NewtonEulerJointR(
            self.joint1_inter.relation())
        self.bar = self.topo.getDynamicalSystem('bar')
        self.post = self.topo.getDynamicalSystem('post')
        self.y = Kernel.SiconosVector(5)
        self.yDoF = Kernel.SiconosVector(1)
        self.jachq = Kernel.SimpleMatrix(1,14)

    def step(self):
        q0 = Kernel.BlockVector(self.bar.q(), self.post.q())
        self.joint1.computeh(0, q0, self.y)
        self.joint1.computehDoF(0, q0, self.yDoF)
        self.joint1.computeJachqDoF(0, self.joint1_inter, q0, self.jachq, 0)
        print('joint angle',self.yDoF)

# Run the simulation
with Hdf5(mode='r+') as io:
    io.run(t0=0,
           T=10,
           h=0.001,
           theta=0.50001,
           Newton_max_iter=1,
           solver=Numerics.SICONOS_FRICTION_3D_NSGS,
           itermax=100000,
           tolerance=1e-10,
           controller=Ctrl(),
           set_external_forces=lambda x: None, # no gravity
           projection_itermax=3,
           projection_tolerance=1e-5,
           projection_tolerance_unilateral=1e-5,
           time_stepping=Kernel.TimeSteppingDirectProjection,
           osi=Kernel.MoreauJeanDirectProjectionOSI,
    )
