"""A few tests for classes and functions from kernel/modelingtools

"""
import numpy as np
import siconos.modeling as mod
import siconos.algebra as alg
 
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
mass[2, 2] = 3.0 / 5 * r * r

weight = np.array([-m * g, 0, 0])
tol = np.finfo(np.double).eps



def test_LagrangianLinearTIDS():
    size = 10
    q0 = alg.SiconosVector(size=size, storage_type=alg.UblasType.dense)
    v0 = alg.SiconosVector(size=size, storage_type=alg.UblasType.dense)
    mass = alg.SimpleMatrix(row=size, col=size, storage_type=alg.UblasType.dense)
    # With numpy array
    ball = mod.LagrangianLinearTIDS(q, v0, mass);
    mg = alg.SiconosVector(size=size)
    mgn = np.array(mg, copy=False)
    mgn[...] = 10.1
    ball.fext = mg
    ball.computeForces(3., q0, v0)
    print(q.data)
    alg.numpy_to_vec(q)
    # assert np.allclose(ball.q(), q, rtol=tol, atol=tol)
    # assert np.allclose(ball.velocity(), v, rtol=tol, atol=tol)
    # assert np.allclose(ball.mass(), mass, rtol=tol, atol=tol)
    # ball.setFExtPtr(weight)
    # assert np.allclose(ball.fExt(), weight, rtol=tol, atol=tol)


def test_NewtonImpactNSL():
    nslaw = mod.NewtonImpactNSL(e=e)
    assert nslaw.e == e
    assert nslaw.size == 1
    nslaw.e = 0.7
    print(nslaw)


def test_LagrangianLinearTIR():
    ndof = 3
    nslawsize = 1
    H = alg.SimpleMatrix(row=nslawsize, col=ndof, storage_type=alg.UblasType.dense)
    b = alg.SiconosVector(size=ndof, storage_type=alg.UblasType.dense)
    #  H = np.array([[1, 0, 0]])
    #     b = np.zeros(1)
    relation = mod.LagrangianLinearTIR(C=H, e=b)
    print(relation)
    #     assert np.allclose(relation.jachq(), H, rtol=tol, atol=tol)


# def test_Nsds():
#     bouncing_ball = K.NonSmoothDynamicalSystem(t0, T)
#     assert bouncing_ball.t0() == t0


if __name__ == '__main__':
    test_LagrangianLinearTIDS()
    test_NewtonImpactNSL()
    test_LagrangianLinearTIR()
