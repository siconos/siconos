# Copyright (C) 2005, 2018 by INRIA
#!/usr/bin/env python

import numpy as np
import os

import siconos.numerics as N
from siconos.tests_setup import working_dir


ztol = 1e-8

def createMLCP_fromFile():
    mlcp=N.MLCP()
    N.mixedLinearComplementarity_newFromFilename(mlcp, os.path.join(working_dir, "data/diodeBridge_mlcp.dat"))
    zsol =np.array([  9.85185185e-01,   9.85185185e-01,  -0.00000000e+00,
         9.85185185e-04,   0.00000000e+00,   0.00000000e+00,
         9.85185185e-04])
    return (mlcp, zsol)

def createMLCP_small():
    # basic interface
    # Murty88, p2
    M = np.array([[2., 1.],
                  [1., 2.]])

    q = np.array([-5., -6.])


    # solution
    zsol = np.array([4. / 3., 7. / 3.])


    # problem
    mlcp=N.MLCP(1,M,q)

    return (mlcp, zsol)

def createMLCP_large():
    zsol =np.array([  9.85185185e-01,   9.85185185e-01,  -0.00000000e+00,
         9.85185185e-04,   0.00000000e+00,   0.00000000e+00,
         9.85185185e-04])


    M = np.array([[  0.00000000e+00,  -1.00000000e-03,   1.00000000e-03,
              0.00000000e+00,   1.00000000e+00,   0.00000000e+00,
              1.00000000e+00],
           [  0.00000000e+00,   1.00000000e-03,  -1.00000000e-03,
             -1.00000000e+00,   0.00000000e+00,  -1.00000000e+00,
              0.00000000e+00],
           [ -1.00250000e+00,   0.00000000e+00,   0.00000000e+00,
              0.00000000e+00,   0.00000000e+00,   1.00000000e+01,
             -1.00000000e+01],
           [  0.00000000e+00,   0.00000000e+00,  -1.00000000e+00,
              0.00000000e+00,   0.00000000e+00,   0.00000000e+00,
              0.00000000e+00],
           [  0.00000000e+00,   1.00000000e+00,   0.00000000e+00,
              0.00000000e+00,   0.00000000e+00,   0.00000000e+00,
              0.00000000e+00],
           [  1.00000000e+00,   0.00000000e+00,  -1.00000000e+00,
              0.00000000e+00,   0.00000000e+00,   0.00000000e+00,
              0.00000000e+00],
           [ -1.00000000e+00,   1.00000000e+00,   0.00000000e+00,
              0.00000000e+00,   0.00000000e+00,   0.00000000e+00,
              0.00000000e+00]])


    q= np.array([[ 0.    ],
           [ 0.    ],
           [ 0.9975],
           [ 0.    ],
           [ 0.    ],
           [ 0.    ],
           [ 0.    ]])

    mlcp = N.MLCP(3, M, q)
    return (mlcp, zsol)


def test_mlcp_enum():
    z = np.array([0., 0.])
    w = np.array([0., 0.])
    mlcp, zsol = createMLCP_small()
    SO=N.SolverOptions(mlcp,N.SICONOS_MLCP_ENUM)
    N.mlcp_driver_init(mlcp, SO)
    info = N.mlcp_enum(mlcp, z, w, SO)
    N.mlcp_driver_reset(mlcp, SO)
    print("z = ", z)
    print("w = ", w)
    assert (np.linalg.norm(z-zsol) <= ztol)
    assert not info


def test_mlcp_enum_large():
    z = np.array([0., 0., 0., 0., 0., 0., 0.])
    w = np.array([0., 0., 0., 0., 0., 0., 0.])
    mlcp, zsol = createMLCP_large()
    SO=N.SolverOptions(mlcp,N.SICONOS_MLCP_ENUM)
    N.mlcp_driver_init(mlcp, SO)
    info = N.mlcp_enum(mlcp, z, w, SO)
    N.mlcp_driver_reset(mlcp, SO)
    print("z = ", z)
    print("w = ", w)
    assert (np.linalg.norm(z-zsol) <= ztol)
    assert not info
#mlcp =0
#N.mixedLinearComplementarity_display(mlcp)

def test_mlcp_enum_large_fromfile():
    z = np.array([0., 0., 0., 0., 0., 0., 0.])
    w = np.array([0., 0., 0., 0., 0., 0., 0.])
    mlcp, zsol = createMLCP_fromFile()
    SO=N.SolverOptions(mlcp,N.SICONOS_MLCP_ENUM)
    N.mlcp_driver_init(mlcp, SO)
    info = N.mlcp_enum(mlcp, z, w, SO)
    N.mlcp_driver_reset(mlcp, SO)
    print("z = ", z)
    print("w = ", w)
    assert (np.linalg.norm(z-zsol) <= ztol)
    assert not info
