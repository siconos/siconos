#!/usr/bin/env python

from Siconos.Kernel import *
from numpy import *

r = 0.1
x = array([1,0,0])
v = array([0,0,0])
mass = eye(3)
mass[2,2]=3./5 * r * r

ds = LagrangianLinearTIDS(x,v,mass)

H = array([[0,0,1]])

nslaw = NewtonImpactNSL(0.9)

relation = LagrangianLinearTIR(H)

inter = Interaction(1, nslaw, relation)
