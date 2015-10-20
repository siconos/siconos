#!/usr/bin/env python

import numpy as np

# import Siconos.Numerics * fails with py.test!
import siconos.numerics as N


NC = 1

M = np.eye(3*NC)

q = np.array([-1., 1., 3.])

mu = np.array([0.1])

z = np.array([0.,0.,0.])

reactions = np.array([0.,0.,0.])

velocities = np.array([0.,0.,0.])


def test_fc3dnsgs():
    N.setNumericsVerbose(2)
    FCP = N.FrictionContactProblem(3,M,q,mu)
    SO=N.SolverOptions(N.SICONOS_FRICTION_3D_NSGS)
    r=N.fc3d_nsgs(FCP, reactions, velocities, SO)
    assert SO.dparam[1] < 1e-10
    assert not r


def test_fc3dlocalac():
    N.setNumericsVerbose(2)
    FCP = N.FrictionContactProblem(3,M,q,mu)
    SO=N.SolverOptions(N.SICONOS_FRICTION_3D_NSN_AC)

    r = N.fc3d_nonsmooth_Newton_AlartCurnier(FCP, reactions, velocities, SO)
    assert SO.dparam[1] < 1e-10
    assert not r


def test_fc3dfischer():
    N.setNumericsVerbose(2)
    FCP = N.FrictionContactProblem(3,M,q,mu)
    SO=N.SolverOptions(N.SICONOS_FRICTION_3D_NSN_FB)

    r = N.fc3d_nonsmooth_Newton_FischerBurmeister(FCP, reactions, velocities, SO)
    assert SO.dparam[1] < 1e-10
    assert not r
