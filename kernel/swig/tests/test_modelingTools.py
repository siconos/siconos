"""A few tests for classes and functions from kernel/modelingtools

"""
#!/usr/bin/env python

import numpy as np
import siconos.kernel as K

t0 = 0
T = 10
r = 0.1
g = 9.81
m = 1
e = 0.9
theta = 0.5
h = 0.005


q = np.array([1, 0, 0])
v = np.array([0, 0, 0])
mass = np.eye(3)
mass[2, 2] = 3. / 5 * r * r

weight = np.array([-m * g, 0, 0])
tol = np.finfo(np.double).eps


def test_LagrangianLinearTIDS():
    ball = K.LagrangianLinearTIDS(q, v, mass)
    assert np.allclose(ball.q(), q, rtol=tol, atol=tol)
    assert np.allclose(ball.velocity(), v, rtol=tol, atol=tol)
    assert np.allclose(ball.mass(), mass, rtol=tol, atol=tol)
    ball.setFExtPtr(weight)
    assert np.allclose(ball.fExt(), weight, rtol=tol, atol=tol)


def test_NewtonImpactNSL():
    nslaw = K.NewtonImpactNSL(e)
    assert nslaw.e() == e


def test_LagrangianLinearTIR():
    H = np.array([[1, 0, 0]])
    b = np.zeros(1)
    relation = K.LagrangianLinearTIR(H, b)
    assert np.allclose(relation.jachq(), H, rtol=tol, atol=tol)


def test_Model():
    bouncing_ball = K.Model(t0, T)
    assert bouncing_ball.t0() == t0

