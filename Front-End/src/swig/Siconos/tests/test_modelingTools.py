#!/usr/bin/env python

from Siconos.Kernel import *
from numpy import *

t0 = 0
T = 10
r = 0.1
g = 9.81
m = 1
e = 0.9
theta = 0.5
h = 0.005


x = array([1,0,0])
v = array([0,0,0])
mass = eye(3)
mass[2,2]=3./5 * r * r

weight = array([-m * g, 0, 0])

ball = LagrangianLinearTIDS(x,v,mass)

ball.setFExtPtr(weight)

H = array([[0,0,1]])

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

#while(s.nextTime() < T):
#    s.computeOneStep()
# ...    
