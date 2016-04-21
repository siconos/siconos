#!/usr/bin/env python

# Siconos-sample, Copyright INRIA 2005-2011.
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
#


import matplotlib
matplotlib.use('Agg')
from matplotlib.pyplot import subplot, title, plot, grid, savefig
from numpy import array, eye, empty, zeros, savetxt
from siconos.kernel import FirstOrderLinearDS, FirstOrderLinearTIR, RelayNSL, \
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
alpha = .01
N = ceil((T-t0)/h)

# matrices
A = zeros((2,2))
x0 = array([10.,10.])
B = 500*array([[alpha,1-alpha],[-(1-alpha),alpha]])
C = eye(2)
D = zeros((2,2))

# dynamical systems
process = FirstOrderLinearDS(x0, A)
myProcessRelation = FirstOrderLinearTIR(C,B)
myProcessRelation.setDPtr(D)

myNslaw = RelayNSL(2)
myNslaw.display()

myProcessInteraction = Interaction(ninter, myNslaw,
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

filippov.initialize(s);

# matrix to save data
dataPlot = empty((N+1,4))
dataPlot[0, 0] = t0
dataPlot[0, 1:3] = process.x()
dataPlot[0, 3] = myProcessInteraction.lambda_(0)[0]

# time loop
k = 1
while(s.hasNextEvent()):
     s.computeOneStep()
     dataPlot[k, 0] = s.nextTime()
     dataPlot[k, 1] = process.x()[0]
     dataPlot[k, 2] = process.x()[1]
     dataPlot[k, 3] = myProcessInteraction.lambda_(0)[0]
     k += 1
     s.nextStep()
     #print s.nextTime()

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
savefig("filipov1.png")

plot(dataPlot[:,1], dataPlot[:,2])
grid()
savefig("filipov2.png")

