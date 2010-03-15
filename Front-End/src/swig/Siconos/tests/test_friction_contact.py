#!/usr/bin/env python

from numpy import *
from Siconos.Numerics import *


NC = 3

M = eye(9)

q = array([-1, 1, 3, -1, 1, 3, -1, 1, 3])

mu = array([0.1, 0.1, 0.1]);

z = array([0,0,0,0,0,0,0,0,0])

setNumericsVerbose(2)

def test_nsgs():
    SO = SolverOptions()
    frictionContact3D_nsgs_setDefaultSolverOptions(SO)
    info = frictionContact3D_driver(FrictionContactProblem(3,M,q,mu),q,z,SO, NumericsOptions())
    print q
    print z

