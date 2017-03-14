"""Tests for python interface of friction contact solvers in numerics.
"""
#!/usr/bin/env python

import numpy as np
import siconos.numerics as sn


NC = 1

M = np.eye(3 * NC)

q = np.array([-1., 1., 3.])

mu = np.array([0.1])

reactions = np.array([0., 0., 0.])
velocities = np.array([0., 0., 0.])
sn.numerics_set_verbose(2)
FCP = sn.FrictionContactProblem(3, M, q, mu)


def solve(problem, solver, options):
    """Solve problem for a given solver
    """
    reactions[...] = 0.0
    velocities[...] = 0.0
    r = solver(problem, reactions, velocities, options)
    assert options.dparam[1] < options.dparam[0]
    assert not r


def test_fc3dnsgs():
    """Non-smooth Gauss Seidel, default
    """
    SO = sn.SolverOptions(sn.SICONOS_FRICTION_3D_NSGS)
    solve(FCP, sn.fc3d_nsgs, SO)


def test_fc3dlocalac():
    """Non-smooth Gauss Seidel, Alart-Curnier as local solver.
    """
    SO = sn.SolverOptions(sn.SICONOS_FRICTION_3D_NSN_AC)
    solve(FCP, sn.fc3d_nonsmooth_Newton_AlartCurnier, SO)


def test_fc3dfischer():
    """Non-smooth Newton, Fischer-Burmeister.
    """
    SO = sn.SolverOptions(sn.SICONOS_FRICTION_3D_NSN_FB)
    solve(FCP, sn.fc3d_nonsmooth_Newton_FischerBurmeister, SO)
