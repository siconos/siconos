#!/usr/bin/env python

# Siconos-sample, Copyright INRIA 2005-2010.
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


from matplotlib.pyplot import *
from Siconos.Kernel import *
from numpy import *

t0 = 0      # start time
T = 10      # end time
h = 0.005   # time step
r = 0.1     # ball radius
g = 9.81    # gravity
m = 1       # ball mass
e = 0.9     # restitution coeficient
theta = 0.5 # theta scheme


x = array([1,0,0])
v = array([0,0,0])
mass = eye(3)
mass[2,2]=3./5 * r * r

weight = array([-m * g, 0, 0])

ball = LagrangianLinearTIDS(x,v,mass)

ball.setFExtPtr(weight)

H = array([[1,0,0]])

nslaw = NewtonImpactNSL(e)

relation = LagrangianLinearTIR(H)

inter = Interaction(1, nslaw, relation)

bouncingBall = Model(t0,T)

bouncingBall.nonSmoothDynamicalSystem().insertDynamicalSystem(ball)

bouncingBall.nonSmoothDynamicalSystem().link(inter,ball);

OSI = Moreau(theta)

OSI.insertDynamicalSystem(ball)

t = TimeDiscretisation(t0,h)

osnspb = LCP()

s = TimeStepping(t)

s.insertIntegrator(OSI)

s.insertNonSmoothProblem(osnspb)

bouncingBall.initialize(s)

N = (T-t0)/h

dataPlot = empty((N+1,5))

lambda_ = inter.lambda_(1)

dataPlot[0, 0] = t0
dataPlot[0, 1] = ball.q()[0]
dataPlot[0, 2] = ball.velocity()[0]
dataPlot[0, 3] = ball.p(2)[0]
dataPlot[0, 4] = inter.lambda_(1)

k = 1
while(s.nextTime() < T):
    s.computeOneStep()

    dataPlot[k, 0] = s.nextTime()
    dataPlot[k, 1] = ball.q()[0]
    dataPlot[k, 2] = ball.velocity()[0]
    dataPlot[k, 3] = ball.p(2)[0]
    dataPlot[k, 4] = inter.lambda_(1)[0]

    k += 1
    s.nextStep()

subplot(411)
title('position')
plot(dataPlot[:,0], dataPlot[:,1])
grid()
subplot(412)
title('velocity')
plot(dataPlot[:,0], dataPlot[:,2])
grid()
subplot(413)
plot(dataPlot[:,0], dataPlot[:,3])
title('reaction')
grid()
subplot(414)
plot(dataPlot[:,0], dataPlot[:,4])
title('lambda')
grid()
show()
