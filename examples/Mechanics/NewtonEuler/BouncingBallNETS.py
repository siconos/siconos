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

from siconos.kernel import NewtonEulerDS, NewtonImpactNSL,\
     NewtonEulerR, NewtonEulerFrom1DLocalFrameR, Interaction, Model,\
     MoreauJeanOSI, TimeDiscretisation, LCP, TimeStepping

from numpy import eye, empty, linalg, savetxt

import math


class BouncingBallR(NewtonEulerFrom1DLocalFrameR):

    def __init__(self, ballRadius):
        self._ballRadius = ballRadius
        NewtonEulerFrom1DLocalFrameR.__init__(self)
        super(BouncingBallR, self).__init__()

    def computeOutput(self, time, interaction, interProp, derivativeNumber):

        #print(interProp.DSlink)
        if derivativeNumber == 0:
            self.computeh(interProp.DSlink[NewtonEulerR.q0], interaction.y(0))
        else:
            super(BouncingBallR, self).computeOutput(t, inter, interProp, derivativeNumber);

    def computeh(self, q, y):
        height = q[0] - self._ballRadius

        y[0] = height

        nnc = [1,0,0]
        self.setnc(nnc)

        ppc1 = [height, q[1], q[2]]
        self.setpc1(ppc1)

        ppc2 = [0.0, q[1], q[2]]
        self.setpc2(ppc2)

t0 = 0      # start time
T = 10.0    # end time
h = 0.005   # time step
r = 0.1     # ball radius
g = 9.81    # gravity
m = 1       # ball mass
e = 0.9     # restitution coeficient
theta = 0.5  # theta scheme

#
# dynamical system
#
x = [1.0, 0, 0, 1.0, 0, 0, 0]  # initial configuration
v = [2.0, 0, 0, 0, 0, 0]  # initial velocity
inertia = eye(3)       # inertia matrix

# the dynamical system
ball = NewtonEulerDS(x, v, m, inertia)

# set external forces
weight = [-m * g, 0, 0]
ball.setFExtPtr(weight)

#
# Interactions
#

# ball-floor
nslaw = NewtonImpactNSL(e)
relation = BouncingBallR(r)

inter = Interaction(nslaw, relation)

#
# Model
#
bouncingBall = Model(t0, T)

# add the dynamical system to the non smooth dynamical system
bouncingBall.nonSmoothDynamicalSystem().insertDynamicalSystem(ball)

# link the interaction and the dynamical system
bouncingBall.nonSmoothDynamicalSystem().link(inter, ball)


#
# Simulation
#

# (1) OneStepIntegrators
OSI = MoreauJeanOSI(theta)

# (2) Time discretisation --
t = TimeDiscretisation(t0, h)

# (3) one step non smooth problem
osnspb = LCP()

# (4) Simulation setup with (1) (2) (3)
s = TimeStepping(t)
s.insertIntegrator(OSI)
s.insertNonSmoothProblem(osnspb)
bouncingBall.setSimulation(s)
# end of model definition

#
# computation
#

# simulation initialization
bouncingBall.initialize()

# the number of time steps
N = (T - t0) / h

# Get the values to be plotted
# ->saved in a matrix dataPlot
dataPlot = empty((N, 16))

#
# numpy pointers on dense Siconos vectors
#
q = ball.q()
v = ball.twist()
p = ball.p(1)
lambda_ = inter.lambda_(1)

#
# initial data
#
dataPlot[0, 0] = t0
dataPlot[0, 1] = q[0]
dataPlot[0, 2] = v[0]
dataPlot[0, 3] = p[0]
dataPlot[0, 4] = lambda_[0]
dataPlot[0, 5] = math.acos(q[3])
dataPlot[0, 6] = linalg.norm(relation.contactForce())
dataPlot[0, 7] = q[0]
dataPlot[0, 8] = q[1]
dataPlot[0, 9] = q[2]
dataPlot[0, 10] = q[3]
dataPlot[0, 11] = q[4]
dataPlot[0, 12] = q[5]
dataPlot[0, 13] = q[6]
dataPlot[0, 14] = v[1]
dataPlot[0, 15] = v[2]
k = 1

# time loop
while(s.hasNextEvent() and k < 2000):
    s.computeOneStep()

    dataPlot[k, 0] = s.nextTime()
    dataPlot[k, 1] = q[0]
    dataPlot[k, 2] = v[0]
    dataPlot[k, 3] = p[0]
    dataPlot[k, 4] = lambda_[0]
    dataPlot[k, 5] = math.acos(q[3])
    dataPlot[k, 6] = linalg.norm(relation.contactForce())
    dataPlot[k, 7] = q[0]
    dataPlot[k, 8] = q[1]
    dataPlot[k, 9] = q[2]
    dataPlot[k, 10] = q[3]
    dataPlot[k, 11] = q[4]
    dataPlot[k, 12] = q[5]
    dataPlot[k, 13] = q[6]
    dataPlot[k, 14] = v[1]
    dataPlot[k, 15] = v[2]
    k = k + 1
    s.nextStep()

savetxt("result-py.dat", dataPlot)
#
# comparison with the reference file
#
from siconos.kernel import SimpleMatrix, getMatrix

ref = getMatrix(SimpleMatrix("resultNETS.ref"))
err = linalg.norm(dataPlot - ref)
print("error w.r.t reference file =", err)

if (err > 1e-12):
    print("Warning. The result is rather different from the reference file.")


#
# plots
#
from matplotlib.pyplot import subplot, title, plot, grid, show

subplot(411)
title('position')
plot(dataPlot[:, 0], dataPlot[:, 1])
grid()
subplot(412)
title('velocity')
plot(dataPlot[:, 0], dataPlot[:, 2])
grid()
subplot(413)
plot(dataPlot[:, 0], dataPlot[:, 3])
title('reaction')
grid()
subplot(414)
plot(dataPlot[:, 0], dataPlot[:, 4])
title('lambda')
grid()
show()
