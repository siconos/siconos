#!/usr/bin/env python3

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

from siconos.kernel import \
    Model, MoreauJeanOSI, TimeDiscretisation, \
    FrictionContact, NewtonImpactFrictionNSL

import siconos.kernel as sk

from siconos.mechanics.proposed import BodyDS, \
    BodyTimeStepping, SiconosBox, SiconosPlane, \
    SiconosContactor, BulletBroadphase, \
    BulletOptions

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
position_init = 3
velocity_init = 0

box = SiconosBox([0,0,0], [1,1,1])
box.setInsideMargin(0.3)
box.setOutsideMargin(0.5)

contactor = SiconosContactor()
contactor.addShape(box)

# A Body Dynamical System : a contactor + a mass + position and velocity
mass = 1.0
body = BodyDS([0, 0, position_init, 1., 0, 0, 0],
              [0, 0, velocity_init, 0., 0., 0.],
              mass)
body.setContactor(contactor)

# set external forces
weight = [0, 0, - mass * g]
body.setFExtPtr(weight)

ground = SiconosContactor()
groundShape = SiconosPlane([0, 0, -.5])
ground.addShape(groundShape)

#
# Model
#
model = Model(t0, T)

# add the dynamical system to the non smooth dynamical system
model.nonSmoothDynamicalSystem().insertDynamicalSystem(body)

#
# Simulation
#

# (1) OneStepIntegrators
osi = MoreauJeanOSI(theta)
osi.insertDynamicalSystem(body)

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
options = BulletOptions()
options.breakingThreshold = 0.5
options.worldScale = 1.0
broadphase = BulletBroadphase(options)

# insert a non smooth law for contactors id 0
# broadphase.insert(nslaw, 0, 0)
# TODO

# (6) Simulation setup with (1) (2) (3) (4) (5)
simulation = BodyTimeStepping(timedisc)
#simulation.setNewtonOptions(1)
simulation.insertIntegrator(osi)
simulation.insertNonSmoothProblem(osnspb)

# simulation initialization
model.initialize(simulation)

broadphase.buildGraph(model)
broadphase.buildGraph(ground)

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

    broadphase.resetStatistics()
    broadphase.updateGraph()
    broadphase.performBroadphase()
    simulation.computeOneStep()

    dataPlot[k, 0] = simulation.nextTime()
    dataPlot[k, 1] = q[2]
    dataPlot[k, 2] = v[2]

    k += 1
    simulation.nextStep()

# #
# # comparison with the reference file
# #
# from siconos.kernel import SimpleMatrix, getMatrix
# from numpy.linalg import norm

# ref = getMatrix(SimpleMatrix("result.ref"))

# print("norm(dataPlot - ref) = {0}".format(norm(dataPlot - ref)))
# if (norm(dataPlot - ref) > 1e-11):
#     print("Warning. The result is rather different from the reference file.")

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
# subplot(313)
# plot(dataPlot[0:k, 0], dataPlot[0:k, 3])
# title('lambda')
# grid()
show()
