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
from siconos.kernel import FirstOrderLinearDS
from siconos.control.simulation import ControlZOHSimulation
from siconos.control.sensor import LinearSensor
from siconos.control.controller import ExplicitLinearSMC

import matplotlib
matplotlib.use('Agg')
from matplotlib.pyplot import subplot, title, plot, grid, savefig, xlabel, ylabel
from numpy import eye, empty, zeros, savetxt
from math import ceil
from numpy.linalg import norm
from matplotlib import rc
import matplotlib.pyplot as plt
import scipy
from scipy import arange

import distutils.spawn
if distutils.spawn.find_executable('latex'):
    rc('text', usetex=True)

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
# Control Simulation
sim = ControlZOHSimulation(t0, T, h)
sim.setSaveOnlyMainSimulation(True)
sim.addDynamicalSystem(processDS)
# Actuator, Sensor & ControlManager
sens = LinearSensor(processDS, sensorC)
sim.addSensor(sens, hControl)
act = ExplicitLinearSMC(sens)
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
savetxt('SMCExampleExplicit-py.dat', dataPlot)

# plot interesting stuff
subplot(211)
ylabel(r'$\sigma$')
xlabel(r't')
plot(dataPlot[:, 0], dataPlot[:, 2])
grid()
subplot(212)
ylabel(r'$\bar{u}^s$')
xlabel(r't')
plt.ylim(-2.1, 2.1)
plot(dataPlot[:, 0], dataPlot[:, 3])
savefig('esmc_sigma_u.png')

subplot(211)
ylabel(r'$\sigma$')
xlabel(r't')
plt.xlim(xmin=.49)
plt.ylim(-0.03, 0.03)
plot(dataPlot[4900:, 0], dataPlot[4900:, 2])
plot([dataPlot[4900, 0], dataPlot[-1, 0]], [0, 0], linewidth=3, color='g', linestyle='dashed')
grid()
subplot(212)
ylabel(r'$\bar{u}^s$')
xlabel(r't')
plt.ylim(-2.1, 2.1)
plt.xlim(xmin=.49)
p1 = plot(dataPlot[4900:, 0], dataPlot[4900:, 3])
#p2 = plot(dataPlot[4900:, 0], np.sin(50*dataPlot[4900:, 0]))
#plt.legend((p1[0], p2[0]), (r'$\bar{u}^s(t)$', r'$-\rho(t)$'), ncol=2)
savefig('esmc_sigma_u_z')

u_z = dataPlot[4900:, 3]
n = len(u_z)
Y = scipy.fft(dataPlot[4900:, 3])/n
k = arange(n)
T = n*h
frq = k/T
frq = frq[list(range(n/2))]
Y = Y[list(range(n/2))]
plot(frq, abs(Y), 'r')
xlabel(r'freq (Hz)')
title(r'Frequency spectrum of $\bar{u}^s$')
savefig('esmc_u_freq.png')
# TODO
# compare with the reference
#ref = getMatrix(SimpleMatrix("result.ref"))
#if (norm(dataPlot - ref[1:,:]) > 1e-12):
#    print("Warning. The result is rather different from the reference file.")
