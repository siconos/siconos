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
import matplotlib.pyplot as plt
from numpy import array, eye, empty, zeros
import numpy as np
from siconos.kernel import FirstOrderLinearDS, RelayNSL, \
NonSmoothDynamicalSystem, Model, TimeDiscretisation, TimeStepping, EulerMoreauOSI, \
Interaction, Relay
from math import ceil


from ZhuravlevTwistingR import ZhuravlevTwistingR


# variables
t0 = 0.0   # start time
T = 10.0     # end time
h = 1.0e-3   # time step
numInter = 2
ninter = 2
theta = 0.5
alpha = .01
N = int((T-t0)/h)

# matrices
A = zeros((2,2))
A[0, 1] = 1

x0 = array([1.,10.])
B = 500*array([[alpha,1-alpha],[-(1-alpha),alpha]])
C = eye(2)
D = zeros((2,2))

# dynamical systems
process = FirstOrderLinearDS(x0, A)
myProcessRelation = ZhuravlevTwistingR(C,B)

myNslaw = RelayNSL(2)
myNslaw.display()

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


#TODO python <- SICONOS_RELAY_LEMKE
# access dparam

osnspb = Relay()
s.insertNonSmoothProblem(osnspb)
s.setComputeResiduY(True)
s.setComputeResiduR(True)

filippov.setSimulation(s)
filippov.initialize()

# matrix to save data
dataPlot = empty((N+1,5))
control = empty((N+1,))
dataPlot[0, 0] = t0
dataPlot[0, 1:3] = process.x()
dataPlot[0, 3] = myProcessInteraction.lambda_(0)[0]
dataPlot[0, 4] = myProcessInteraction.lambda_(0)[1]
# time loop
k = 1
while(s.hasNextEvent()):
     s.newtonSolve(1e-14, 30)
     dataPlot[k, 0] = s.nextTime()
     dataPlot[k, 1] = process.x()[0]
     dataPlot[k, 2] = process.x()[1]
     dataPlot[k, 3] = myProcessInteraction.lambda_(0)[0]
     dataPlot[k, 4] = myProcessInteraction.lambda_(0)[1]
     control[k] = process.r()[1]
     k += 1
     s.nextStep()
     #print s.nextTime()

# save to disk
np.savetxt('ZI_Twisting.txt', dataPlot)
np.savetxt('ZI_Twisting_u.txt', control)
# plot interesting stuff
plt.subplot(411)
plt.title('s')
plt.plot(dataPlot[:,0], dataPlot[:,1])
plt.grid()
plt.subplot(412)
plt.title('v')
plt.plot(dataPlot[:,0], dataPlot[:,2])
plt.grid()
plt.subplot(413)
plt.plot(dataPlot[:,0], dataPlot[:,3])
plt.title('lambda1')
plt.grid()
plt.subplot(414)
plt.plot(dataPlot[:,0], dataPlot[:,4])
plt.title('lambda2')
plt.grid()
plt.savefig('ZI_Twisting_all.png')

plt.plot(dataPlot[:,1], dataPlot[:,2])
plt.xlabel('s')
plt.xlabel('v')
plt.grid()
plt.savefig('ZI_Twisting_sv.png')

plt.plot(dataPlot[:,3], dataPlot[:,4])
plt.xlabel('lambda1')
plt.xlabel('lambda2')
plt.grid()
plt.savefig('ZI_Twisting_lambdas.png')

pos = np.abs(dataPlot[:,1])
velocity = (1-myProcessRelation._kappa*np.sign(dataPlot[:,1]*dataPlot[:,2]))*dataPlot[:, 2]*np.sign(dataPlot[:,1])

plt.subplot(311)
plt.title('position')
plt.plot(dataPlot[:,0], pos)
plt.grid()
plt.subplot(312)
plt.title('velocity')
plt.plot(dataPlot[:,0], velocity)
plt.grid()
plt.subplot(313)
plt.title('control input')
plt.plot(dataPlot[:,0], control)
plt.grid()
plt.savefig('ZI_Twisting_pv.png')

indx = np.nonzero(dataPlot[:, 0]>30)
ttt = dataPlot[indx, 0].flatten()

plt.subplot(311)
plt.title('position')
plt.plot(ttt, pos[indx])
plt.grid()
plt.subplot(312)
plt.title('velocity')
plt.plot(ttt, velocity[indx])
plt.grid()
plt.subplot(313)
plt.title('control input')
plt.plot(ttt, control[indx])
plt.grid()
plt.savefig('ZI_Twisting_pv_z.png')
