#!/usr/bin/env python

from siconos.tests_setup import working_dir
import os
import siconos.kernel as sk
import numpy as np

ndof = 3
nn = 2 * ndof
x0 = np.zeros(nn, dtype=np.float64)
q0 = np.zeros(ndof, dtype=np.float64)
v0 = np.zeros(ndof, dtype=np.float64)
q_ne = np.zeros(nn + 1, dtype=np.float64)
inertia = np.eye(ndof, dtype=np.float64)
a_mat = np.asarray(np.random.random((nn, nn)), dtype=np.float64)

# List of all existing DS classes and standard initialisation paramaters
ds_classes = {
    sk.FirstOrderNonLinearDS: (x0,),
    sk.FirstOrderLinearDS: (x0, a_mat),
    sk.FirstOrderLinearTIDS: (x0, a_mat),
    sk.LagrangianDS: (q0, v0, inertia),
    sk.LagrangianLinearTIDS: (q0, v0, inertia),
    sk.NewtonEulerDS: (q_ne, x0, 1., inertia)}


def test_ds_interface():
    """Tests methods that should be available
    for all dynamical systems
    """
    time = 1.

    for args in ds_classes:
        class_name = args
        attr = ds_classes[args]
        ds = class_name(*attr)
        print(class_name, attr[0])
        ds.initializeNonSmoothInput(1)
        ds.resetAllNonSmoothParts()
        ds.resetNonSmoothPart(1)
        ds.initRhs(time)
        ds.computeRhs(time)
        assert ds.number() >= 0
        assert ds.rhs() is not None
        ds.update(time)
        ds.resetToInitialState()
        ds.display()
        ds.computeJacobianRhsx(time)
        assert ds.jacobianRhsx() is not None


def test_first_order_nlds():
    """Build and test first order non linear ds
    """
    ndof = 3
    x0 = np.zeros(ndof, dtype=np.float64)
    ds = sk.FirstOrderNonLinearDS(x0)
    ds.display()
    time = 1.2
    ds.computeRhs(time)
    ds.computeJacobianRhsx(time)
    assert ds.dimension() == ndof
    assert not ds.isLinear()
    assert np.allclose(ds.x0(), x0)
    assert np.allclose(ds.x(), x0)
    assert np.allclose(ds.rhs(), 0.)
    ds.computef(time, ds.x())
    assert ds.f() is None
    ds.initRhs(time)
    assert ds.jacobianfx() is None
    assert np.allclose(ds.jacobianRhsx(), 0.)


def test_first_order_lds():
    """Build and test first order linear
    and time-invariant coeff. ds
    """
    ndof = 3
    x0 = np.random.random(ndof)
    time = 1.2
    ds_list = []
    a_mat = np.random.random((ndof, ndof))
    b_vec = np.random.random((ndof, ))
    ds_list.append(sk.FirstOrderLinearDS(x0))
    ds_list.append(sk.FirstOrderLinearDS(x0, a_mat))
    ds_list.append(sk.FirstOrderLinearDS(x0, a_mat, b_vec))

    for ds in ds_list:
        assert ds.isLinear()
        assert ds.dimension() == ndof
        assert np.allclose(ds.x0(), x0)
        assert np.allclose(ds.x(), x0)
        assert np.allclose(ds.r(), 0.)

        rhs = np.zeros_like(ds.x())
        jac_ref = np.zeros((ndof, ndof), dtype=np.float64)
        if isinstance(ds.A(), np.ndarray):
            jac_ref += a_mat
            rhs += np.dot(a_mat, ds.x())
        if isinstance(ds.b(), np.ndarray):
            rhs += ds.b()
        ds.computef(time, ds.x())
        if ds.f() is not None:
            assert np.allclose(rhs, ds.f())

        ds.initRhs(time)
        assert np.allclose(rhs, ds.rhs())
        if ds.A() is not None:
            assert np.allclose(ds.jacobianRhsx(), jac_ref)
            assert np.allclose(ds.jacobianRhsx(), ds.jacobianfx())


def test_first_order_ltids():
    """Build and test first order linear
    and time-invariant coeff. ds
    """
    time = 1.2
    ds_list = []
    b_vec = np.random.random((nn, ))
    ds_list.append(sk.FirstOrderLinearTIDS(x0, a_mat))
    ds_list.append(sk.FirstOrderLinearTIDS(x0, a_mat, b_vec))

    for ds in ds_list:
        assert ds.isLinear()
        assert ds.dimension() == nn
        assert np.allclose(ds.x0(), x0)
        assert np.allclose(ds.x(), x0)
        assert np.allclose(ds.r(), 0.)

        rhs = np.dot(a_mat, ds.x())
        if ds.b() is not None:
            rhs += ds.b()
        ds.initRhs(time)
        assert np.allclose(rhs, ds.rhs())
        assert np.allclose(ds.jacobianRhsx(), a_mat)
        assert np.allclose(ds.jacobianRhsx(), ds.jacobianfx())


def test_lagrangian_ds():
    """Build and test lagrangian ds
    """
    q0[...] = [1, 2, 3]
    v0[...] = [4, 5, 6]

    mass = np.asarray(np.diag([1, 2, 3]), dtype=np.float64)
    ds = sk.LagrangianDS(q0, v0, mass)
    ec = ds.computeKineticEnergy()
    assert ec == 87.
    assert ds.dimension() == ndof
    assert np.allclose(ds.mass(), mass)


def test_lagrangian_tids():
    """Build and test lagrangian linear and time-invariant ds
    """
    q0[...] = [1, 2, 3]
    v0[...] = [4, 5, 6]

    mass = np.asarray(np.diag([1, 2, 3]), dtype=np.float64)
    stiffness = np.zeros((ndof, ndof), dtype=np.float64)
    stiffness.flat[...] = np.arange(9)
    damping = np.zeros_like(stiffness)
    damping.flat[...] = np.arange(9, 18)
    ds = sk.LagrangianLinearTIDS(q0, v0, mass, stiffness, damping)
    ec = ds.computeKineticEnergy()
    assert ec == 87.
    assert ds.dimension() == ndof
    assert np.allclose(ds.mass(), mass)
    assert np.allclose(ds.K(), stiffness)
    assert np.allclose(ds.C(), damping)
    q = ds.q()
    v = ds.velocity()
    fref = -np.dot(stiffness, q)
    fref -= np.dot(damping, v)
    time = 0.3
    ds.computeForces(time, q, v)
    assert np.allclose(fref, ds.forces())
    ds.computeJacobianqForces(time)
    assert np.allclose(stiffness, ds.jacobianqForces())
    ds.computeJacobianqDotForces(time)
    assert np.allclose(damping, ds.jacobianqDotForces())
