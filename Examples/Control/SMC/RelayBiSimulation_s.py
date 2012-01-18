#!/usr/bin/env python

# Siconos-sample, Copyright INRIA 2005-2012.
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

from Siconos.Kernel import *
from matplotlib.pyplot import subplot, title, plot, grid, show
from numpy import array, eye, empty, zeros, savetxt
from math import ceil
from numpy.linalg import norm

# variable declaration
ndof = 2   # Number of degrees of freedom of your system
t0 = 0.0   # start time
T = 1    # end time
h = 1.0e-4  # time step for simulation
hControl = 1.0e-2 # time step for control
Xinit = 1.0 # initial position
theta = 0.5
N = ceil((T-t0)/h + 10) # number of time steps
outputSize = 3 # number of variable to store at each time step

# Matrix declaration
A = zeros((ndof,ndof))
x0 = [Xinit, -Xinit]
sensorC = eye(ndof)
sensorD = zeros((ndof,ndof))
Csurface = eye(ndof)

# Simple check
if h > hControl:
    print "hControl must be bigger than h"
    exit(1)

# Declaration of the Dynamical System
controlProcess = ControlFirstOrderLinearDS(t0, T, h, x0, A)
processDS = controlProcess.processDS();
processDS.setComputebFunction("RelayPlugin.so","computeB")
controlProcess.initialize()
controlModel = controlProcess.model()
processSimulation = controlProcess.simulation()

# time discretisation
tSensor = TimeDiscretisation(t0, hControl)
tActuator = TimeDiscretisation(t0, hControl)
# Actuator, Sensor & ControlManager
control = ControlManager(controlModel)
sens = linearSensor(100, tSensor, controlModel, sensorC, sensorD)
control.addSensorPtr(sens)
act = linearSMC(101, tActuator, controlModel)
act.addSensorPtr(sens)
control.addActuatorPtr(act)

# Initialization
control.initialize()
act.setCsurfacePtr(Csurface)
# This is not working right now
#eventsManager = s.eventsManager()

# Matrix for data storage
dataPlot = empty((3*(N+1), outputSize))
dataPlot[0, 0] = t0
dataPlot[0, 1] = processDS.x()[0]
dataPlot[0, 2] = processDS.x()[1]

# Main loop
k = 1
while(processSimulation.nextTime() < T):
    processSimulation.computeOneStep()
    dataPlot[k, 0] = processSimulation.nextTime()
    dataPlot[k, 1] = processDS.x()[0]
    dataPlot[k, 2] = processDS.x()[1]
    k += 1
    processSimulation.nextStep()
#    print processSimulation.nextTime()
# Resize
dataPlot.resize(k, outputSize)
# Save to disk
savetxt('RelayBiSimulation-py.dat', dataPlot)
# Plot interesting data
subplot(211)
title('x1')
plot(dataPlot[:,0], dataPlot[:,1])
grid()
subplot(212)
title('x2')
plot(dataPlot[:,0], dataPlot[:,2])
grid()
show()
# TODO
# compare with the reference
#ref = getMatrix(SimpleMatrix("result.ref"))
#if (norm(dataPlot - ref[1:,:]) > 1e-12):
#    print("Warning. The result is rather different from the reference file.")
