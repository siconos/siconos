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
from siconos.kernel import FirstOrderLinearDS, Model, TimeDiscretisation, \
    TimeStepping, ZeroOrderHoldOSI
from siconos.control.simulation import ControlManager
from siconos.control.sensor import LinearSensor
from siconos.control.controller import LinearSMCOT2
from matplotlib.pyplot import subplot, title, plot, grid, show, xlabel, ylabel
import matplotlib.pyplot as plt
from numpy import eye, empty, zeros, savetxt
from math import ceil

# variable declaration
ndof = 2   # Number of degrees of freedom of your system
t0 = 0.0   # start time
T = 1    # end time
h = 1.0e-4  # time step for simulation
hControl = 1.0e-2  # time step for control
Xinit = 1.0  # initial position
theta = 0.5
N = 2*ceil((T-t0)/h)  # number of time steps
outputSize = 5  # number of variable to store at each time step

# Matrix declaration
A = zeros((ndof, ndof))
x0 = [Xinit, -Xinit]
sensorC = eye(ndof)
Csurface = [[0, 1.0]]

# Simple check
if h > hControl:
    print "hControl must be bigger than h"
    exit(1)

# Declaration of the Dynamical System
processDS = FirstOrderLinearDS(x0, A)
processDS.setComputebFunction("RelayPlugin", "computeB")
# Model
process = Model(t0, T)
process.nonSmoothDynamicalSystem().insertDynamicalSystem(processDS)
# time discretisation
processTD = TimeDiscretisation(t0, h)
tSensor = TimeDiscretisation(t0, hControl)
tActuator = TimeDiscretisation(t0, hControl)
# Creation of the Simulation
processSimulation = TimeStepping(processTD, 0)
processSimulation.setName("plant simulation")
# Declaration of the integrator
processIntegrator = ZeroOrderHoldOSI(processDS)
processSimulation.insertIntegrator(processIntegrator)
# Actuator, Sensor & ControlManager
control = ControlManager(process)
sens = LinearSensor(tSensor, processDS, sensorC)
control.addSensorPtr(sens)
act = LinearSMCOT2(tActuator, processDS)
act.setCsurfacePtr(Csurface)
act.addSensorPtr(sens)
control.addActuatorPtr(act)

# Initialization
process.initialize(processSimulation)
control.initialize()
# This is not working right now
#eventsManager = s.eventsManager()

# Matrix for data storage
dataPlot = empty((N+1, outputSize))
#dataPlot[0, 0] = processDS.t0()
dataPlot[0, 0] = t0
dataPlot[0, 1] = processDS.x()[0]
dataPlot[0, 2] = processDS.x()[1]
dataPlot[0, 3] = processDS.z()[0]
dataPlot[0, 4] = processDS.z()[1]

# Main loop
k = 1
while(processSimulation.hasNextEvent()):
    processSimulation.computeOneStep()
    dataPlot[k, 0] = processSimulation.nextTime()
    dataPlot[k, 1] = processDS.x()[0]
    dataPlot[k, 2] = processDS.x()[1]
    dataPlot[k, 3] = processDS.z()[0]
    dataPlot[k, 4] = processDS.z()[1]
    k += 1
    print processSimulation.nextTime()
    processSimulation.nextStep()
# Resize matrix
dataPlot.resize(k, outputSize)
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
title('u1')
plot(dataPlot[:, 0], dataPlot[:, 3])
grid()
subplot(414)
title('u2')
plot(dataPlot[:, 0], dataPlot[:, 4])
show()

subplot(211)
p1 = plot(dataPlot[:, 0], dataPlot[:, 2])
plt.ylabel('$\sigma$')
plt.xlabel('t')
grid()
subplot(212)
p2 = plot(dataPlot[:, 0], dataPlot[:, 4])
plt.ylabel('u')
plt.xlabel('t')
show()

# TODO
# compare with the reference
#ref = getMatrix(SimpleMatrix("result.ref"))
#if (norm(dataPlot - ref[1:,:]) > 1e-12):
#    print("Warning. The result is rather different from the reference file.")
