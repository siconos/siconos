#!/usr/bin/env python

# Siconos is a program dedicated to modeling, simulation and control
# of non smooth dynamical systems.
#
# Copyright 2018 INRIA.
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

import sys
import matplotlib
matplotlib.use('Agg')
from matplotlib.pyplot import subplot, title, plot, grid, savefig, figure
from numpy import array, eye, empty, zeros, savetxt, loadtxt, linalg
from siconos.kernel import FirstOrderLinearDS, FirstOrderLinearTIR, RelayNSL,\
    NonSmoothDynamicalSystem, TimeDiscretisation, TimeStepping, EulerMoreauOSI, \
Interaction, Relay
from math import ceil


import MyR, MyNonLinearR


# variables
t0 = 0.0   # start time
T = 0.1     # end time
h = 1.0e-3   # time step
numInter = 2
ninter = 2
theta = 0.5
alpha = .01
N = int((T-t0)/h)

# matrices
A = zeros((2,2))
x0 = array([3.,3.])
B = 500*array([[alpha,1-alpha],[-(1-alpha),alpha]])
C = eye(2)
D = zeros((2,2))

# dynamical systems
process = FirstOrderLinearDS(x0, A)
#myProcessRelation = MyR.MyR(C,B)
myProcessRelation = MyNonLinearR.MyNonLinearR(C,B)
#myProcessRelation.setDPtr(D)

myNslaw = RelayNSL(2)
#myNslaw.display()

myProcessInteraction = Interaction(myNslaw,
        myProcessRelation)




filippov = NonSmoothDynamicalSystem(t0,T)
filippov.insertDynamicalSystem(process)
filippov.link(myProcessInteraction,process)


#myProcessRelation.computeJachx(0, x0, x0 , x0, C)

td = TimeDiscretisation(t0, h)
s = TimeStepping(filippov, td)


myIntegrator = EulerMoreauOSI(theta)

s.insertIntegrator(myIntegrator)


#TODO python <- SICONOS_RELAY_LEMKE
# access dparam

osnspb = Relay()
s.insertNonSmoothProblem(osnspb)
s.setComputeResiduY(True)
s.setComputeResiduR(True)


# matrix to save data
dataPlot = empty((N+1,4))
dataPlot[0, 0] = t0
dataPlot[0, 1:3] = process.x()
dataPlot[0, 3] = myProcessInteraction.lambda_(0)[0]

# time loop
k = 1
while(s.hasNextEvent()):
     if (k%10):
          sys.stdout.write('.')

     s.computeOneStep()
     #osnspb.display()
     dataPlot[k, 0] = s.nextTime()
     dataPlot[k, 1] = process.x()[0]
     dataPlot[k, 2] = process.x()[1]
     dataPlot[k, 3] = myProcessInteraction.lambda_(0)[0]
     k += 1
     s.nextStep()
     #print s.nextTime()
sys.stdout.write('\n')
# save to disk
savetxt('Filippov_NL_py.dat', dataPlot)

dataRef = loadtxt('Filippov_NL_py.dat')

print('Comparison with reference file -  error = ',linalg.norm(dataPlot-dataRef) )
assert(linalg.norm(dataPlot-dataRef) <= 1e-12)


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
savefig("Filipov_NL1.png")

figure()
plot(dataPlot[:,1], dataPlot[:,2])
grid()
savefig("Filipov_NL2.png")
