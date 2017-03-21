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
import math

t0 = 0.0
T = 100   # Total simulation time
h_step = 1.0e-2  # Time step
xinit = 3.0*math.sqrt(2.0)/(2.0*math.pi)    # initial voltage
Modeltitle = "RelayOscillator"

withPlot=True
if withPlot:
    import matplotlib
    matplotlib.use('Agg')
    from matplotlib.pyplot import subplot, title, plot, grid, savefig

from siconos.kernel import FirstOrderLinearDS, FirstOrderLinearTIR, \
                           RelayNSL, Interaction,\
                           Model, EulerMoreauOSI, TimeDiscretisation, Relay,  \
                           TimeStepping,NonSmoothDynamicalSystem

#
# dynamical system
#
init_state = [0.0,xinit,0]

A = [[0,          1.0,     0.0],
     [0.0,        0.0,     1.0],
     [0.0,       -3.0,    -2.0]]

LSRelayOscillator=FirstOrderLinearDS(init_state, A)

#
# Interactions
#

C = [[1.0, 0.0,   0.0]]

D = [[0.0 ]]

B = [[0.0],
     [0.0],
     [1.0]]

LTIRRelayOscillator=FirstOrderLinearTIR(C,B)
LTIRRelayOscillator.setDPtr(D)

nslaw=RelayNSL(1)
InterRelayOscillator=Interaction(nslaw,LTIRRelayOscillator)


#
# Model
#
RelayOscillator=Model(t0,T,Modeltitle)

#   add the dynamical system in the non smooth dynamical system
myNSDS = NonSmoothDynamicalSystem()
myNSDS.insertDynamicalSystem(LSRelayOscillator)


#   link the interaction and the dynamical system
myNSDS.link(InterRelayOscillator,LSRelayOscillator)

RelayOscillator.setNonSmoothDynamicalSystemPtr(myNSDS)


#
# Simulation
#

# (1) OneStepIntegrators
theta = 0.5
aOSI = EulerMoreauOSI(theta)
# (2) Time discretisation
aTiDisc = TimeDiscretisation(t0,h_step)

# (3) Non smooth problem
aRelay = Relay()

# (4) Simulation setup with (1) (2) (3)
aTS = TimeStepping(aTiDisc,aOSI,aRelay)

# end of model definition

#
# computation
#

# simulation initialization
RelayOscillator.setSimulation(aTS)
RelayOscillator.initialize()

k = 0
h = aTS.timeStep();
print("Timestep : ",h)
# Number of time steps
N = (int)((T-t0)/h)
print("Number of steps : ",N)

# Get the values to be plotted
# ->saved in a matrix dataPlot

from numpy import empty
dataPlot = empty([N+1,8])

x = LSRelayOscillator.x()
print("Initial state : ",x)
y = InterRelayOscillator.y(0)
print("First y : ",y)
lambda_ = InterRelayOscillator.lambda_(0)

while (k < N):
    aTS.computeOneStep()
    #aLCP.display()
    dataPlot[k, 0] = aTS.nextTime()
    #  inductor voltage
    dataPlot[k, 1] = x[0]
    dataPlot[k, 2] = x[1]
    dataPlot[k, 3] = x[2]
    dataPlot[k, 4] = y[0]
    dataPlot[k, 5] = lambda_[0]


    k += 1
    if k%1000==0:
        print("step =", k, " < ", N)
    aTS.nextStep()

if (withPlot) :
    #
    # plots
    #
    subplot(511)
    title('x1')
    plot(dataPlot[0:k-1,0], dataPlot[0:k-1,1])
    grid()
    subplot(512)
    title('x2')
    plot(dataPlot[0:k-1,0], dataPlot[0:k-1,2])
    grid()
    subplot(513)
    title('x3')
    plot(dataPlot[0:k-1,0], dataPlot[0:k-1,3])
    subplot(514)
    title('y')
    plot(dataPlot[0:k-1,0], dataPlot[0:k-1,4])
    grid()
    subplot(515)
    title('lambda')
    plot(dataPlot[0:k-1,0], dataPlot[0:k-1,5])
    grid()
    savefig("relay_oscillator.png")

