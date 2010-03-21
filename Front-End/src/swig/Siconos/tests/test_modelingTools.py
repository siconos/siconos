#!/usr/bin/env python

from matplotlib.pyplot import *
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


q = array([1,0,0])
v = array([0,0,0])
mass = eye(3)
mass[2,2]=3./5 * r * r

weight = array([-m * g, 0, 0])

def test_LagrangianLinearTIDS():
    ball = LagrangianLinearTIDS(q,v,mass)
    assert (ball.q().all() == q.all())
    assert (ball.velocity().all() == v.all())
    assert (ball.mass().all() == mass.all())

    ball.setFExtPtr(weight)


def test_NewtonImpactNSL():
    nslaw = NewtonImpactNSL(e)
    assert(nslaw.e() == e)

def test_LagrangianLinearTIR():
    H = array([[1,0,0]])
    relation = LagrangianLinearTIR(H)
    assert(relation.H() == H)

def test_Model():
    bouncingBall = Model(t0,T)
    assert (bouncingBall.t0() == t0)


