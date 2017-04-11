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

import siconos.kernel as sk
from siconos.control.simulation import ControlManager
from siconos.control.sensor import LinearSensor
from siconos.control.controller import PID

import matplotlib
matplotlib.use('Agg')
from matplotlib.pyplot import subplot, title, plot, grid, savefig
from numpy import array, eye, empty, zeros, savetxt
from math import ceil
from numpy.linalg import norm

# variable declaration
t0 = 0.0   # start time
T = 100.0  # end time
h = 0.05  # time step
xFinal = 0.0  # target position
theta = 0.5
N = int((T - t0) / h + 10)  # number of time steps
outputSize = 5  # number of variable to store at each time step

# Matrix declaration
A = zeros((2, 2))
A[0, 1] = 1
B = [[0], [1]]
x0 = [10., 10.]
C = [[1., 0]]  # we have to specify ndmin=2, so it's understood as
K = [.25, .125, 2]

# Declaration of the Dynamical System
doubleIntegrator = sk.FirstOrderLinearTIDS(x0, A)
# Model
process = sk.Model(t0, T)
process.nonSmoothDynamicalSystem().insertDynamicalSystem(doubleIntegrator)
# Declaration of the integrator
OSI = sk.EulerMoreauOSI(theta)
# time discretisation
t = sk.TimeDiscretisation(t0, h)
tSensor = sk.TimeDiscretisation(t0, h)
tActuator = sk.TimeDiscretisation(t0, h)
s = sk.TimeStepping(t, 0)
s.insertIntegrator(OSI)

# Actuator, Sensor & ControlManager
control = ControlManager(s)
sens = LinearSensor(doubleIntegrator, C)
control.addSensorPtr(sens, tSensor)
act = PID(sens)
act.setB(B)
control.addActuatorPtr(act, tActuator)

# Initialization
process.setSimulation(s)
process.initialize()
control.initialize(process)
act.setRef(xFinal)
act.setK(K)
act.setDeltaT(h)
# This is not working right now
#eventsManager = s.eventsManager()

# Matrix for data storage
dataPlot = zeros((3 * N + 3, outputSize))
dataPlot[0, 0] = process.t0()
dataPlot[0, 1] = doubleIntegrator.x()[0]
dataPlot[0, 2] = doubleIntegrator.x()[1]
if doubleIntegrator.b() is not None:
    dataPlot[0, 3] = doubleIntegrator.b()[0]
    dataPlot[0, 4] = doubleIntegrator.b()[1]

# Main loop
k = 1
while(s.hasNextEvent()):
    if (s.eventsManager().nextEvent() == 1):
        s.computeOneStep()
        dataPlot[k, 0] = s.nextTime()
        dataPlot[k, 1] = doubleIntegrator.x()[0]
        dataPlot[k, 2] = doubleIntegrator.x()[1]
        if doubleIntegrator.b() is not None:
            dataPlot[k, 3] = doubleIntegrator.b()[0]
            dataPlot[k, 4] = doubleIntegrator.b()[1]
        k += 1
    s.nextStep()
# Save to disk
savetxt('output.txt', dataPlot)
# Plot interesting data
subplot(211)
title('position')
plot(dataPlot[:, 0], dataPlot[:, 1])
grid()
subplot(212)
title('velocity')
plot(dataPlot[:, 0], dataPlot[:, 2])
grid()
savefig("pid.png")
# TODO
# compare with the reference
#ref = getMatrix(SimpleMatrix("result.ref"))
#if (norm(dataPlot - ref[1:,:]) > 1e-12):
#    print("Warning. The result is rather different from the reference file.")
