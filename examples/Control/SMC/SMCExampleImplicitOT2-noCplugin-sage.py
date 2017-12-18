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

# We have first to import sage, so it won't shadow some function (like plot...)
# You need the env. variable DOT_SAGE to point to the right location !


import sys

try:
    from sage.all import *
    # this is needed since sage uses mpfr by default ...
    RealNumber = float
    Integer = int
except ImportError:
    print('sage is not installed, exiting')
    sys.exit(0)

# Other import
from siconos.kernel import FirstOrderLinearDS, getMatrix
from siconos.control.simulation import ControlZOHSimulation
from siconos.control.sensor import LinearSensor
from siconos.control.controller import LinearSMCOT2

import matplotlib
matplotlib.use('Agg')
from matplotlib.pyplot import subplot, title, plot, grid, savefig
from numpy import array, eye, empty, zeros, savetxt
from math import ceil, sin
from numpy.linalg import norm


# Some stupid symbolic computations
x = var('x')
g = vector((-cos(50*x), cos(50*x)))/50
f = g.diff(x)

# Derive our own version of FirstOrderLinearDS

class MyFOLDS(FirstOrderLinearDS):
    def computeb(self, time):
        tmpz = self.z()
        # XXX fix this !
        if len(tmpz) != 2:
            print("DEBUG z has length ", len(tmpz))
            return
        # XXX we need to find a smarter way to do things here
        # we need to convert from vector (sage) to arrayish
        u = array(f(x=time).list(), dtype = float) + tmpz
        self.setbPtr(u)

# variable declaration
ndof = 2   # Number of degrees of freedom of your system
t0 = 0.0   # start time
T = 1    # end time
h = 1.0e-4  # time step for simulation
hControl = 1.0e-2 # time step for control
Xinit = 1.0 # initial position
theta = 0.5
N = ceil((T-t0)/h + 10) # number of time steps
outputSize = 5 # number of variable to store at each time step

# Matrix declaration
A = zeros((ndof, ndof))
x0 = [Xinit, -Xinit]
sensorC = eye(ndof)
Csurface = [[0, 1.0]]

# Simple check
if h > hControl:
    print("hControl must be bigger than h")
    exit(1)

# Declaration of the Dynamical System
processDS = MyFOLDS(x0, A)
# XXX b is not automatically created ...
processDS.setbPtr([0, 0])

# Control Simulation
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
savetxt('SMCExampleImplicitOT2-noCplugin-sage-py.dat', dataPlot)
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
savefig('ismcOT2_x_u.png')

# compare with the reference
ref = getMatrix(SimpleMatrix("SMCExampleImplicitOT2-py.ref"))
print("%19e" % norm(dataPlot - ref))
if (norm(dataPlot - ref) > 1e-12):
    print(dataPlot - ref)
    print("Warning. The result is rather different from the reference file.")
