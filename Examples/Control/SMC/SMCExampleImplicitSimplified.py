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

from siconos.kernel import ControlFirstOrderLinearDS, LinearSMC, getMatrix, \
    TimeDiscretisation, SimpleMatrix
from siconos.control.sensor import LinearSensor
from matplotlib.pyplot import subplot, title, plot, grid, show
from numpy import eye, zeros, savetxt
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
A = zeros((ndof, ndof))
x0 = [Xinit, -Xinit]
sensorC = eye(ndof)
Csurface = eye(ndof)
Brel = eye(ndof)
Drel = eye(ndof)

# Simple check
if h > hControl:
    print "hControl must be bigger than h"
    exit(1)

# Declaration of the Dynamical System
controlProcess = ControlFirstOrderLinearDS(t0, T, h, x0, A)
processDS = controlProcess.processDS()
processDS.setComputebFunction("RelayPlugin", "computeB")
controlProcess.initialize()

# time discretisation
tSensor = TimeDiscretisation(t0, hControl)
tActuator = TimeDiscretisation(t0, hControl)
# Actuator, Sensor & ControlManager
control = controlProcess.CM()
sens = LinearSensor(tSensor, processDS, sensorC)
control.addSensorPtr(sens)
act = LinearSMC(tActuator, processDS)
act.setCsurfacePtr(Csurface)
act.setBPtr(Brel)
act.setDPtr(Drel)
act.addSensorPtr(sens)
control.addActuatorPtr(act)

# Initialization
control.initialize()

# Run the simulation
controlProcess.run()
# get the results
tmpData = controlProcess.data()
dataPlot = tmpData[:, 0:3]
# Save to disk
savetxt('SMCExampleImplicitSimplified-py.dat', dataPlot)
# Plot interesting data
subplot(211)
title('x1')
plot(dataPlot[:, 0], dataPlot[:, 1])
grid()
subplot(212)
title('x2')
plot(dataPlot[:, 0], dataPlot[:, 2])
grid()
show()
# compare with the reference
ref = getMatrix(SimpleMatrix("SMCExampleImplicit.ref"))
print('\n')
print(norm(dataPlot - ref))
if (norm(dataPlot - ref) > 1e-12):
    print("Warning. The result is rather different from the reference file.")
