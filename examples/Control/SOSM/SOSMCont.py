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


import matplotlib
matplotlib.use('Agg')
from matplotlib.pyplot import subplot, title, plot, grid, savefig
from numpy import eye, empty, zeros, savetxt
from siconos.kernel import FirstOrderNonLinearDS, FirstOrderLinearTIR, RelayNSL, \
NonSmoothDynamicalSystem, Model, TimeDiscretisation, TimeStepping, EulerMoreauOSI, \
Interaction, Relay
from math import ceil

# variables
t0 = 0.0   # start time
T = 10      # end time
h = 1.0e-3   # time step
numInter = 2
ninter = 2
theta = 0.5
N = (int)((T-t0)/h)
mu1 = 2
mu2 = 3

# matrices
A = zeros((2,2))
A[0,1] = 1
#x0 = array([10.,0.])
x0 = ([10.,0.])
#x0 = 10

#B = 500*array([[0,0],[mu2,mu1]])
B = ([0,0],[mu2,mu1])
#B = mu2

C = eye(2)
#C = 1
D = zeros((2,2))
#D = 0

# dynamical systems
process = FirstOrderNonLinearDS(x0)
process.setComputeFFunction('PluginF', 'computef1')
process.setComputeJacobianfxFunction('PluginF', 'computeJacf1')


#process = FirstOrderNonLinearDS(x0,'PluginF:computef1','PluginF:computeJacf1'  )

process.display()

myProcessRelation = FirstOrderLinearTIR(C,B)
myProcessRelation.setDPtr = D

myNslaw = RelayNSL(2)
myNslaw.display()

nameInter = 'processInteraction'
myProcessInteraction = Interaction(myNslaw,
        myProcessRelation)
myNSDS = NonSmoothDynamicalSystem()
myNSDS.insertDynamicalSystem(process)
myNSDS.link(myProcessInteraction,process)


filippov = Model(t0,T)
filippov.setNonSmoothDynamicalSystemPtr(myNSDS)

td = TimeDiscretisation(t0, h)
s = TimeStepping(td)

myIntegrator = EulerMoreauOSI(theta)
s.insertIntegrator(myIntegrator)

print('initialization')
#TODO python <- SICONOS_RELAY_LEMKE
# access dparam

osnspb = Relay()
s.insertNonSmoothProblem(osnspb)
filippov.setSimulation(s)
filippov.initialize()

print('end of initialization')
# matrix to save data
dataPlot = empty((N+1,4))
dataPlot[0, 0] = t0
dataPlot[0, 1:3] = process.x()
dataPlot[0, 3] = myProcessInteraction.lambda_(0)[0]

# time loop
k = 1
while(s.hasNextEvent()):
     #print 'iteration k'
     s.computeOneStep()
     dataPlot[k, 0] = s.nextTime()
     dataPlot[k, 1] = process.x()[0]
     dataPlot[k, 2] = process.x()[1]
     dataPlot[k, 3] = myProcessInteraction.lambda_(0)[0]
     k += 1
     s.nextStep()

# save to disk
savetxt('output.txt', dataPlot)
# plot interesting stuff
subplot(311)
title('position')
plot(dataPlot[:,0], dataPlot[:,1])
grid()
subplot(312)
title('velocity')
plot(dataPlot[:,0], dataPlot[:,2])
grid()
subplot(313)
plot(dataPlot[:,0], dataPlot[:,3])
title('lambda')
grid()
savefig("SOSM1.png")

plot(dataPlot[:,1], dataPlot[:,2])
grid()
savefig("SOSM1.png")
