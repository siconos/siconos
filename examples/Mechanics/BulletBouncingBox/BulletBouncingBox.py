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

do_plot = True
try:
    import matplotlib
except:
    do_plot = False
if do_plot:
    import os, sys
    if sys.platform=='linux' and (not 'DISPLAY' in os.environ
                                  or len(os.environ['DISPLAY'])==0):
        matplotlib.use('Agg')
    from matplotlib.pyplot import \
        subplot, title, plot, grid, show, savefig, ylim

from siconos.kernel import \
    Model, MoreauJeanOSI, TimeDiscretisation, \
    FrictionContact, NewtonImpactFrictionNSL, TimeStepping

import siconos.kernel as sk

from siconos.mechanics.collision.bullet import \
     SiconosBulletCollisionManager

from siconos.mechanics.collision import \
    SiconosBox, SiconosPlane, BodyDS, SiconosContactor, SiconosContactorSet

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

# a box shape
box1 = SiconosBox(1.0, 1.0, 1.0)

# A Bullet Dynamical System : a shape + a mass (1.0) + position and velocity
body = BodyDS([0, 0, position_init, 1., 0, 0, 0],
              [0, 0, velocity_init, 0., 0., 0.],
              1.0)

# Add the shape, wrapped in a SiconosContactor, to the body's
# contactor set.
body.contactors().push_back(SiconosContactor(box1))

# set external forces
weight = [0, 0, -body.scalarMass() * g]
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
osi = MoreauJeanOSI(theta)

ground = SiconosPlane()
groundOffset = [0,0,-0.5,1,0,0,0]

# (2) Time discretisation --
timedisc = TimeDiscretisation(t0, h)

# (3) one step non smooth problem
osnspb = FrictionContact(3)

osnspb.numericsSolverOptions().iparam[0] = 1000
osnspb.numericsSolverOptions().dparam[0] = 1e-5
osnspb.setMaxSize(16384)
osnspb.setMStorageType(1)
osnspb.setNumericsVerboseMode(False)

# keep previous solution
osnspb.setKeepLambdaAndYState(True)


# (4) non smooth law
nslaw = NewtonImpactFrictionNSL(0.8, 0., 0., 3)

# (5) broadphase contact detection
broadphase = SiconosBulletCollisionManager()

# insert a non smooth law for contactors id 0
broadphase.insertNonSmoothLaw(nslaw, 0, 0)

# The ground is a static object
# we give it a group contactor id : 0
scs = SiconosContactorSet()
scs.append(SiconosContactor(ground))
broadphase.insertStaticContactorSet(scs, groundOffset)

# (6) Simulation setup with (1) (2) (3) (4) (5)
simulation = TimeStepping(timedisc)
simulation.insertInteractionManager(broadphase)

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

    simulation.computeOneStep()

    dataPlot[k, 0] = simulation.nextTime()
    dataPlot[k, 1] = q[2]
    dataPlot[k, 2] = v[2]

    #if (broadphase.collisionWorld().getDispatcher().getNumManifolds() > 0):
    if (broadphase.statistics().new_interactions_created +
        broadphase.statistics().existing_interactions_processed) > 0:
        if bouncingBox.nonSmoothDynamicalSystem().topology().\
          numberOfIndexSet() == 2:
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

ref = getMatrix(SimpleMatrix("result.ref"))

print("norm(dataPlot - ref) = {0}".format(norm(dataPlot - ref)))
if (norm(dataPlot - ref) > 1e-11):
    print("Warning. The result is rather different from the reference file.")


#
# plots
#

if do_plot:
    subplot(511)
    title('position')
    plot(dataPlot[0:k, 0], dataPlot[0:k, 1])
    y = ylim()
    plot(ref[0:k, 0], ref[0:k, 1])
    ylim(y)
    grid()
    subplot(513)
    title('velocity')
    plot(dataPlot[0:k, 0], dataPlot[0:k, 2])
    y = ylim()
    plot(ref[0:k, 0], ref[0:k, 2])
    ylim(y)
    grid()
    subplot(515)
    plot(dataPlot[0:k, 0], dataPlot[0:k, 3])
    y = ylim()
    plot(ref[0:k, 0], ref[0:k, 3])
    ylim(y)
    title('lambda')
    grid()
    savefig('result.png')
    show()
