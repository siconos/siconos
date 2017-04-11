#!/usr/bin/env python
# -*- coding: utf-8 -*-
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
#-----------------------------------------------------------------------
#
#  DiodeBridgeCapFilter  : sample of an electrical circuit involving :
#  - a 1st linear dynamical system LSDiodeBridge1 consisting of
#        an LC oscillator (1 muF , 10 mH)
#  - a non smooth system : a 4 diodes bridge used as a full wave rectifier
#        of the supplied voltage across the LC oscillator, providing power
#    to the resistor load of the 2nd dynamical system
#      - a 2nd linear dynamical system LSDiodeBridge2 consisting of
#        a filtering capacitor in parallel with a load resistor
#
#  Expected behavior :
#  The initial state (Vc = 10 V , IL = 0) of the oscillator provides
#      an initial energy.
#  The oscillator period is 2 Pi sqrt(LC) ~ 0,628 ms.
#      The non smooth system is a full wave rectifier :
#  each phase (positive and negative) of the oscillation allows current
#      to flow in a constant direction through the load.
#      The capacitor filter acts as a tank providing energy to the load resistor
#      when the voltage across the oscillator weakens.
#      The load resistor consumes energy : the oscillation damps.
#
#  State variables LSDiodeBridge1:
#  - the voltage across the capacitor (or inductor)
#  - the current through the inductor
#
#  State variable LSDiodeBridge2:
#  - the voltage across the filtering capacitor
#
#  The interaction between the two dynamical systems is defined by :
#  - complementarity laws between diodes current and voltage. Depending on
#        the diode position in the bridge, y stands for the reverse voltage across
#    the diode or for the diode current (see figure in the template file)
#  - a linear time invariant relation between the state variables and y and
#    lambda (derived from the Kirchhoff laws)
#
#-----------------------------------------------------------------------

t0 = 0.0
T = 5.0e-3       # Total simulation time
h_step = 1.0e-6  # Time step
Lvalue = 1e-2    # inductance
Cvalue = 1e-6    # capacitance
Rvalue = 1e3     # resistance
Cfilt  = 300.0e-9 # filtering capacitor
VinitLS1 = 10.0   # initial voltage LC oscillator
VinitLS2 = 0.0    # initial voltage Cfilt

Modeltitle = "DiodeBridge"

withPlot = True
if (withPlot):
    import matplotlib
    matplotlib.use('Agg')
    from matplotlib.pyplot import subplot, title, plot, grid, savefig

from siconos.kernel import FirstOrderLinearDS, FirstOrderLinearTIR, \
                           ComplementarityConditionNSL, Interaction,\
                           Model, EulerMoreauOSI, TimeDiscretisation, LCP,  \
                           TimeStepping

#
# dynamical system
#
init_stateLS1 = [VinitLS1, 0]

LS1_A = [[0,          -1.0/Cvalue],
     [1.0/Lvalue, 0          ]]

LS1DiodeBridgeCapFilter = FirstOrderLinearDS(init_stateLS1, LS1_A)

init_stateLS2 = [VinitLS2]

LS2_A = [[-1.0/(Rvalue*Cfilt)]]

LS2DiodeBridgeCapFilter = FirstOrderLinearDS(init_stateLS2, LS2_A)

#
# Interactions
#

C = [[0.,   0., 1.0],
     [0,    0., 0.0],
     [-1.,  0., 1.0],
     [1.,   0., 0.0]]

D = [[0.0, -1.0, 0.,  0.],
     [1.0,  0.0,  1., -1.],
     [0.0,       -1.,         0.,  0.],
     [0.,         1.,         0.,  0.]]

B = [[0.,        0., -1./Cvalue, 1./Cvalue],
     [0.,        0.,  0.,        0.       ],
     [1.0/Cfilt,        0.,  1.0/Cfilt,        0.       ]]

LTIRDiodeBridgeCapFilter = FirstOrderLinearTIR(C, B)
LTIRDiodeBridgeCapFilter.setDPtr(D)

nslaw = ComplementarityConditionNSL(4)
InterDiodeBridgeCapFilter = Interaction(nslaw, LTIRDiodeBridgeCapFilter)

#
# Model
#
DiodeBridgeCapFilter = Model(t0, T, Modeltitle)

#   add the dynamical system in the non smooth dynamical system
DiodeBridgeCapFilter.nonSmoothDynamicalSystem().insertDynamicalSystem(LS1DiodeBridgeCapFilter)
DiodeBridgeCapFilter.nonSmoothDynamicalSystem().insertDynamicalSystem(LS2DiodeBridgeCapFilter)

#   link the interaction and the dynamical system
DiodeBridgeCapFilter.nonSmoothDynamicalSystem().link(InterDiodeBridgeCapFilter, LS1DiodeBridgeCapFilter, LS2DiodeBridgeCapFilter)

#
# Simulation
#

# (1) OneStepIntegrators
theta = 0.5
gamma = 1.0
aOSI = EulerMoreauOSI(theta, gamma)
aOSI.setUseGammaForRelation(True)

# (2) Time discretisation
aTiDisc = TimeDiscretisation(t0, h_step)

# (3) Non smooth problem
aLCP = LCP()

# (4) Simulation setup with (1) (2) (3)
aTS = TimeStepping(aTiDisc, aOSI, aLCP)

# end of model definition

#
# computation
#

# simulation initialization
DiodeBridgeCapFilter.setSimulation(aTS)
DiodeBridgeCapFilter.initialize()

k = 0
h = aTS.timeStep()
print("Timestep : ", h)
# Number of time steps
N = int((T - t0) / h)
print("Number of steps : ", N)

# Get the values to be plotted
# ->saved in a matrix dataPlot

from numpy import zeros
dataPlot = zeros([N-1, 10])

x = LS1DiodeBridgeCapFilter.x()
print("Initial state : ", x)
y = InterDiodeBridgeCapFilter.y(0)
print("First y : ", y)
lambda_ = InterDiodeBridgeCapFilter.lambda_(0)

# For the initial time step:
# time

#  inductor voltage
dataPlot[k, 1] = x[0]

# inductor current
dataPlot[k, 2] = x[1]

# diode R1 current
dataPlot[k, 3] = lambda_[0]

# diode R1 voltage
dataPlot[k, 4] = - y[0]

# diode F2 voltage
dataPlot[k, 5] = - lambda_[1]

# diode F1 current
dataPlot[k, 6] = lambda_[2]





k += 1
while (k < N-1):
    aTS.computeOneStep()
    #aLCP.display()
    dataPlot[k, 0] = aTS.nextTime()
    #  inductor voltage
    dataPlot[k, 1] = x[0]
    # inductor current
    dataPlot[k, 2] = x[1]
    # diode R1 current
    dataPlot[k, 3] =    lambda_[0]
    # diode R1 voltage
    dataPlot[k, 4] = - y[0]
    # diode F2 voltage
    dataPlot[k, 5] = - lambda_[1]
    # diode F1 current
    dataPlot[k, 6] = lambda_[2]
    k += 1
    aTS.nextStep()

# comparison with reference file
from siconos.kernel import SimpleMatrix, getMatrix
from numpy.linalg import norm

ref = getMatrix(SimpleMatrix("DiodeBridgeCapFilter.ref"))

error = norm(dataPlot[:,0:6] - ref[:,0:6])
print("error = " , error)

#assert (error < 1e-09)
withRef = False
if (withPlot):
    #
    # plots
    #
    subplot(411)
    title('inductor voltage')
    plot(dataPlot[0:k - 1, 0], dataPlot[0:k - 1, 1])
    if (withRef):
        plot(ref[0:k - 1, 0], ref[0:k - 1, 1])
    grid()
    subplot(412)
    title('inductor current')
    plot(dataPlot[0:k - 1, 0], dataPlot[0:k - 1, 2])
    if (withRef):
        plot(ref[0:k - 1, 0], ref[0:k - 1, 2])
    grid()
    subplot(413)
    title('diode R1 (blue) and F2 (green) voltage')
    plot(dataPlot[0:k - 1, 0], -dataPlot[0:k - 1, 4])
    plot(dataPlot[0:k - 1, 0], dataPlot[0:k - 1, 5])
    if (withRef):
        plot(ref[0:k - 1, 0], -ref[0:k - 1, 4])
        plot(ref[0:k - 1, 0], ref[0:k - 1, 5])
    grid()

    subplot(414)
    title('resistor voltage')
    plot(dataPlot[0:k - 1, 0], -dataPlot[0:k - 1, 4] - dataPlot[0:k - 1, 5]  )
    if (withRef):
        plot(dataPlot[0:k - 1, 0], -ref[0:k - 1, 4] - ref[0:k - 1, 5]  )
    grid()
    savefig("diode_bridge_capfilter_tgs.png")
