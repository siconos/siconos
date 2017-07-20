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
from siconos.kernel import FirstOrderLinearDS, getMatrix, SimpleMatrix
from siconos.control.simulation import ControlZOHSimulation
from siconos.control.sensor import LinearSensor
from siconos.control.controller import LinearSMCOT2

import matplotlib
matplotlib.use('Agg')
from matplotlib.pyplot import subplot, title, plot, grid, savefig
import matplotlib.pyplot as plt
from numpy import eye, zeros, savetxt
from numpy.linalg import norm
from math import ceil

# variable declaration
ndof = 2   # Number of degrees of freedom of your system
t0 = 0.0   # start time
T = 1    # end time
h = 1.0e-4  # time step for simulation
hControl = 1.0e-2  # time step for control
Xinit = 1.0  # initial position
theta = 0.5
N = 2*int(ceil((T-t0)/h))  # number of time steps
outputSize = 5  # number of variable to store at each time step

# Matrix declaration
A = zeros((ndof, ndof))
x0 = [Xinit, -Xinit]
sensorC = eye(ndof)
Csurface = [[0, 1.0]]
Brel = [[0], [2]]

# Simple check
if h > hControl:
    print("hControl must be bigger than h")
    exit(1)

# Declaration of the Dynamical System
processDS = FirstOrderLinearDS(x0, A)
processDS.setComputebFunction("RelayPlugin", "computeB")
# Control simulation
sim = ControlZOHSimulation(t0, T, h)
sim.setSaveOnlyMainSimulation(True)
sim.addDynamicalSystem(processDS)
# Actuator, Sensor & ControlManager
sens = LinearSensor(processDS, sensorC)
sim.addSensor(sens, hControl)
act = LinearSMCOT2(sens)
act.setCsurface(Csurface)
act.setB(Brel)
sim.addActuator(act, hControl)

# Initialization
sim.initialize()

# Run simulation
sim.run()

# Get data
dataPlot = sim.data()

# Save to disk
savetxt('SMCExampleImplicitOT2-py.dat', dataPlot)
# Plot interesting data
subplot(411)
title('x1')
plot(dataPlot[:, 0], dataPlot[:, 1])
grid()
subplot(412)
title('x2')
plot(dataPlot[:, 0], dataPlot[:, 2])
grid()
subplot(413)
title('u')
plot(dataPlot[:, 0], dataPlot[:, 3])
grid()
savefig('ismcOT2_x_u')

subplot(211)
p1 = plot(dataPlot[:, 0], dataPlot[:, 2])
plt.ylabel('$\sigma$')
plt.xlabel('t')
grid()
subplot(212)
p2 = plot(dataPlot[:, 0], dataPlot[:, 3])
plt.ylabel('u')
plt.xlabel('t')
savefig('ismcOT2_sigma_u')

# TODO
# compare with the reference
ref = getMatrix(SimpleMatrix("SMCExampleImplicitOT2-py.ref"))
if (norm(dataPlot - ref) > 1e-12):
    print("Warning. The result is rather different from the reference file.")
