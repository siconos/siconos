#!/usr/bin/env python

# Siconos is a program dedicated to modeling, simulation and control
# of non smooth dynamical systems.
#
# Copyright 2016 INRIA.
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
# http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.
#
#

from siconos.kernel import \
    Model, MoreauJeanOSI, TimeDiscretisation, \
    FrictionContact, NewtonImpactFrictionNSL

import siconos.kernel as sk

from siconos.mechanics.collision.bullet import \
     btConvexHullShape, btVector3, btCollisionObject, \
     btBoxShape, btMatrix3x3, \
     BulletSpaceFilter, \
     BulletWeightedShape, BulletDS, BulletTimeStepping

from numpy import zeros
from numpy.linalg import norm

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

def makeBox(pos=position_init, vel=velocity_init):
    if (True):
        box = btConvexHullShape()
        box.addPoint(btVector3(-1.0, 1.0, -1.0))
        box.addPoint(btVector3(-1.0, -1.0, -1.0))
        box.addPoint(btVector3(-1.0, -1.0, 1.0))
        box.addPoint(btVector3(-1.0, 1.0, 1.0))
        box.addPoint(btVector3(1.0, 1.0, 1.0))
        box.addPoint(btVector3(1.0, 1.0, -1.0))
        box.addPoint(btVector3(1.0, -1.0, -1.0))
        box.addPoint(btVector3(1.0, -1.0, 1.0))
    else:
        box = btBoxShape(btVector3(1.0, 1.0, 1.0))

    # a bullet shape with a mass (1.0)
    box1 = BulletWeightedShape(box, 1.0)

    # A Bullet Dynamical System : a shape + a mass + position and velocity
    body = BulletDS(box1,
                    [0, 0, pos, 1., 0, 0, 0],
                    [0, 0, vel, 0., 0., 0.])

    # set external forces
    weight = [0, 0, - box1.mass() * g]
    body.setFExtPtr(weight)

    return body

# Initial box body
body = makeBox()

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
osi = MoreauJeanOSI(theta)

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

osnspb.numericsSolverOptions().iparam[0] = 1000
osnspb.numericsSolverOptions().dparam[0] = 1e-5
osnspb.setMStorageType(1)
osnspb.setNumericsVerboseMode(False)

# keep previous solution
osnspb.setKeepLambdaAndYState(True)


# (4) non smooth law
nslaw = NewtonImpactFrictionNSL(0.8, 0., 0., 3)

# (5) broadphase contact detection
broadphase = BulletSpaceFilter(bouncingBox)

# insert a non smooth law for contactors id 0
broadphase.insert(nslaw, 0, 0)

# add multipoint iterations to gather at least 3 contact points and
# avoid object penetration
broadphase.collisionConfiguration().setConvexConvexMultipointIterations()
broadphase.collisionConfiguration().setPlaneConvexMultipointIterations()

# The ground is a static object
# we give it a group contactor id : 0
broadphase.addStaticObject(ground, 0)

# (6) Simulation setup with (1) (2) (3) (4) (5)
simulation = BulletTimeStepping(timedisc)
#simulation.setNewtonOptions(1)
simulation.insertIntegrator(osi)
simulation.insertNonSmoothProblem(osnspb)


# simulation initialization
bouncingBox.setSimulation(simulation)
bouncingBox.initialize()

# Get the values to be plotted
# ->saved in a matrix dataPlot

N = int((T - t0) / h)
dataPlot = zeros((N+1, 4))

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

    # Add a second box dynamically to the simulation
    if k == 100:
        broadphase.addDynamicObject(makeBox(pos=3), simulation)

    broadphase.buildInteractions(bouncingBox.currentTime())
    simulation.computeOneStep()

    dataPlot[k, 0] = simulation.nextTime()
    dataPlot[k, 1] = q[2]
    dataPlot[k, 2] = v[2]

    if (broadphase.collisionWorld().getDispatcher().getNumManifolds() > 0):
        if bouncingBox.nonSmoothDynamicalSystem().topology().\
          numberOfIndexSet() > 1:
            index1 = sk.interactions(simulation.indexSet(1))
            if (len(index1) == 4):
                dataPlot[k, 3] = norm(index1[0].lambda_(1)) + \
                norm(index1[1].lambda_(1)) + norm(index1[2].lambda_(1)) + \
                norm(index1[3].lambda_(1))

    k += 1
    simulation.nextStep()

#
# comparison with the reference file
#
from siconos.kernel import SimpleMatrix, getMatrix
from numpy.linalg import norm

ref = getMatrix(SimpleMatrix("result_dynamic.ref"))

print("norm(dataPlot - ref) = {0}".format(norm(dataPlot - ref)))
if (norm(dataPlot - ref) > 1e-11):
    print("Warning. The result is rather different from the reference file.")


#
# plots
#

from matplotlib.pyplot import subplot, title, plot, grid, show

subplot(511)
title('position')
plot(dataPlot[0:k, 0], dataPlot[0:k, 1])
grid()
subplot(513)
title('velocity')
plot(dataPlot[0:k, 0], dataPlot[0:k, 2])
grid()
subplot(515)
plot(dataPlot[0:k, 0], dataPlot[0:k, 3])
title('lambda')
grid()
show()
