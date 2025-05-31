"""A few tests for classes and functions from kernel/modelingtools

"""

import numpy as np
import siconos.modeling as sm


def test_lagrangianDS():
    print("start test_lagrangianDS")
    ndof = 3
    q0 = np.array([1, 0, 0], dtype=np.float64)
    v0 = np.zeros_like(q0)
    ball = sm.LagrangianDS(q0, v0)

    fext = np.zeros_like(q0)
    fext[:] = 122
    ball.setConstantFext(fext)
    assert np.allclose(ball.fext()[:], fext[:])
    fext[1] = 12
    assert np.allclose(ball.fext()[:], fext[:])

    fext = 453
    assert ball.fext()[0] == ball.fext()[2] == 122
    assert ball.fext()[1] == 12
    print("end test_lagrangianDS")

def fint_func(vel, pos, time, result):
    result[0] = vel[1] * np.cos(time)
    result[:] = vel * time + pos

def mass_func(pos, result):
    for i in range(pos.size):
        result[:,i] = i * pos

def test_lagrangianDS_compute_mass():
    print("start test_lagrangianDS_compute_mass")
    ndof = 3
    q0 = np.zeros(ndof, dtype=np.float64)
    v0 = np.zeros_like(q0)
    q0[:] = np.arange(1, ndof + 1)
    v0[:] = np.arange(4, ndof + 4)
    ball = sm.LagrangianDS(q0, v0)
    ball.setComputeMassFunction(mass_func)

    q = ball.q()
    v = ball.velocity()
    p = ball.p(1)
    time = 1.0

    q[:] = 1.
    v[:] = 2.
    p[:] = 3.

    ref = np.zeros((ndof, ndof),dtype=np.float64)
    ref[:]=1.
    mass_func(q,ref)

    ball.computeMass(q)
    assert np.allclose(ball.mass, ref)

def test_lagrangianDS_compute():
    print("start test_lagrangianDS_compute")

    ndof = 3
    q0 = np.zeros(ndof, dtype=np.float64)
    v0 = np.zeros_like(q0)
    q0[:] = np.arange(1, ndof + 1)
    v0[:] = np.arange(4, ndof + 4)
    ball = sm.LagrangianDS(q0, v0)
    ball.setComputeFintFunction(fint_func)

    q = ball.q()
    v = ball.velocity()
    p = ball.p(1)
    time = 1.0

    q[:] = 1.
    v[:] = 2.
    p[:] = 3.
    ref = np.zeros(q.size)
    fint_func(v, q, time, ref)

    ball.computeFint(v, q, 1.0)
    assert np.allclose(ball.fint, ref)


if __name__ == "__main__":
    test_lagrangianDS()
    test_lagrangianDS_compute()
    test_lagrangianDS_compute_mass()
