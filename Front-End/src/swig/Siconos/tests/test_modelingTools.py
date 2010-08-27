#!/usr/bin/env python


from numpy import *
import Siconos.Kernel as K

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

def equalv(v1,v2):
    return (linalg.norm(v1-v2) <= finfo(double).eps)
    

def test_LagrangianLinearTIDS():
    ball = K.LagrangianLinearTIDS(q,v,mass)
    assert (equalv(ball.q(),q))
    assert (equalv(ball.velocity(),v))
    assert (equalv(ball.mass(),mass))

    ball.setFExtPtr(weight)

    assert(equalv(ball.fExt(),weight))


def test_NewtonImpactNSL():
    nslaw = K.NewtonImpactNSL(e)
    assert(nslaw.e() == e)

def test_LagrangianLinearTIR():
    H = array([[1,0,0]])
    relation = K.LagrangianLinearTIR(H)
    assert(equalv(relation.jachq(),H))

def test_Model():
    bouncingBall = K.Model(t0,T)
    assert (bouncingBall.t0() == t0)

def test_display():
    ball = K.LagrangianLinearTIDS(q,v,mass)
    ball.display()

def test_number():
    ball = K.LagrangianLinearTIDS(q,v,mass)
    print ball.number()

