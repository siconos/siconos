#!/usr/bin/env python

# Siconos-sample, Copyright INRIA 2005-2013.
# Siconos is a program dedicated to modeling, simulation and control
# of non smooth dynamical systems.
# Siconos is a free software; you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation; either version 2 of the License, or
# (at your option) any later version.
# Siconos is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with Siconos; if not, write to the Free Software
# Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301  USA
#
# Contact: Vincent ACARY, siconos-team@lists.gforge.fr
#

from Siconos.Kernel import \
     Model, Moreau, TimeDiscretisation,\
     FrictionContact, NewtonImpactFrictionNSL

from Siconos.Mechanics.contactDetection import \
    btConvexHullShape, btVector3, btCollisionObject, \
    btBoxShape, btMatrix3x3, \
    BulletSpaceFilter, \
    BulletWeightedShape, BulletDS, BulletTimeStepping

from numpy import empty

t0 = 0       # start time
T = 20       # end time
h = 0.005    # time step

g = 9.81     # gravity

theta = 0.5  # theta scheme

#
# dynamical system
#
position_init = 10
velocity_init = 0

box = btConvexHullShape()
box.addPoint(btVector3(-1.0, 1.0, -1.0))
box.addPoint(btVector3(-1.0, 1.0, -1.0))
box.addPoint(btVector3(-1.0, -1.0, -1.0))
box.addPoint(btVector3(-1.0, -1.0, 1.0))
box.addPoint(btVector3(-1.0, 1.0, 1.0))
box.addPoint(btVector3(1.0, 1.0, 1.0))
box.addPoint(btVector3(1.0, 1.0, -1.0))
box.addPoint(btVector3(1.0, -1.0, -1.0))
box.addPoint(btVector3(1.0, -1.0, 1.0))

# a bullet shape with a mass (1.0)
box1 = BulletWeightedShape(box, 1.0)

# A Bullet Dynamical System : a shape + a mass + position and velocity
body = BulletDS(box1,
                [0, 0, position_init, 1., 0, 0, 0],
                [0, 0, velocity_init, 0., 0., 0.])

# set external forces
weight = [0, 0, - box1.mass() * g ]
body.setFExtPtr(weight)

#
# Model
#
bouncingBox = Model(t0, T)

# add the dynamical system to the non smooth dynamical system
bouncingBox.nonSmoothDynamicalSystem().insertDynamicalSystem(body)

#
# Simulation
#

# (1) OneStepIntegrators
osi = Moreau(theta)
osi.insertDynamicalSystem(body)


ground = btCollisionObject()
ground.setCollisionFlags(btCollisionObject.CF_STATIC_OBJECT)
groundShape = btBoxShape(btVector3(30, 30, .5))
basis = btMatrix3x3()
basis.setIdentity()
ground.getWorldTransform().setBasis(basis)
ground.setCollisionShape(groundShape)
ground.getWorldTransform().getOrigin().setZ(-.5)

# (2) Time discretisation --
timedisc = TimeDiscretisation(t0, h)

# (3) one step non smooth problem
osnspb = FrictionContact(3)

# keep previous solution
osnspb.setKeepLambdaAndYState(True)


# (4) non smooth law
nslaw = NewtonImpactFrictionNSL(0.8, 0., 0., 3)

# (5) broadphase contact detection
aabbmax = btVector3(100, 100, 100)
aabbmin = btVector3(-100, -100, -100)
broadphase = BulletSpaceFilter(bouncingBox, nslaw, aabbmin, aabbmax)

broadphase.collisionWorld().addCollisionObject(ground)
broadphase.addStaticObject(ground)
broadphase.addStaticShape(groundShape)


# (6) Simulation setup with (1) (2) (3) (4) (5)
simulation = BulletTimeStepping(timedisc, broadphase)
simulation.insertIntegrator(osi)
simulation.insertNonSmoothProblem(osnspb)


# simulation initialization
bouncingBox.initialize(simulation)

# Get the values to be plotted
# ->saved in a matrix dataPlot

N = (T - t0) / h
dataPlot = empty((N, 7))

#
# numpy pointers on dense Siconos vectors
#
q = body.q()
v = body.velocity()

#
# initial data
#
dataPlot[0, 0] = t0
dataPlot[0, 1] = q[2]
dataPlot[0, 2] = v[2]

k = 1

# time loop
while(simulation.hasNextEvent()):

    broadphase.buildInteractions(bouncingBox.currentTime())

    simulation.computeOneStep()

    dataPlot[k, 0] = simulation.nextTime()
    dataPlot[k, 1] = q[2]
    dataPlot[k, 2] = v[2]

    if (broadphase.collisionWorld().getDispatcher().getNumManifolds() > 0):
        index1 = simulation.indexSet(1)
        if (index1.size() > 0):
            # output reaction forces here
            # to be done
            pass

    k += 1
    simulation.nextStep()

#
# comparison with the reference file
#
from Siconos.Kernel import SimpleMatrix, getMatrix
from numpy.linalg import norm

ref = getMatrix(SimpleMatrix("result.ref"))

# to be done...
if (norm(dataPlot - ref) > 1e-12):
    print("Warning. The result is rather different from the reference file.")


#
# plots
#

from matplotlib.pyplot import subplot, title, plot, grid, show

subplot(211)
title('position')
plot(dataPlot[0:k, 0], dataPlot[0:k, 1])
grid()
subplot(212)
title('velocity')
plot(dataPlot[0:k, 0], dataPlot[0:k, 2])
grid()
#subplot(414)
#plot(dataPlot[:, 0], dataPlot[:, 4])
#title('lambda')
#grid()
show()
