#!/usr/bin/env python

from numpy import *

# import Siconos.Numerics * fails with py.test!
import Siconos.Numerics as N


NC = 3

M = eye(9)

q = array([-1, 1, 3, -1, 1, 3, -1, 1, 3])

mu = array([0.1, 0.1, 0.1]);

z = array([0,0,0,0,0,0,0,0,0])

reactions = array([0,0,0,0,0,0,0,0,0])

velocities = array([0,0,0,0,0,0,0,0,0])


def test_fc3dnsgs():
    N.setNumericsVerbose(2)
    FCP = N.FrictionContactProblem(3,M,q,mu)
    SO=N.SolverOptions()
    N.frictionContact3D_setDefaultSolverOptions(SO,N.SICONOS_FRICTION_3D_NSGS)
    r=N.frictionContact3D_nsgs(FCP, reactions, velocities, SO)
    assert not r


