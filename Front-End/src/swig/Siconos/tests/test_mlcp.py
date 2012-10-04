#!/usr/bin/env python

from numpy import *

# import Siconos.Numerics * fails with py.test!
import Siconos.Numerics as N

# basic interface
# Murty88, p2
M = array([[2., 1.],
           [1., 2.]])

q = array([-5., -6.])

z = array([0., 0.])

w = array([0., 0.])

# solution
zsol = array([4./3., 7./3.])
wsol = array([0. , 0.])

# problem
mlcp=N.MLCP(1,M,q)

ztol = 1e-8

def test_mlcp_enum():
    SO=N.SolverOptions(mlcp,N.SICONOS_MLCP_ENUM)
    N.mlcp_driver_init(mlcp, SO)
    info = N.mlcp_enum(mlcp, z, w, SO)
    N.mlcp_driver_reset(mlcp, SO)
    print "z = ", z
    print "w = ", w
    assert (linalg.norm(z-zsol) <= ztol)
    assert not info

z = array([0., 0., 0., 0., 0., 0.,0.])
zsol =array([  9.85185185e-01,   9.85185185e-01,  -0.00000000e+00,
         9.85185185e-04,   0.00000000e+00,   0.00000000e+00,
         9.85185185e-04]) 

w = array([0., 0., 0., 0., 0., 0.,0.])

M = array([[  0.00000000e+00,  -1.00000000e-03,   1.00000000e-03,
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


q= array([[ 0.    ],
       [ 0.    ],
       [ 0.9975],
       [ 0.    ],
       [ 0.    ],
       [ 0.    ],
       [ 0.    ]])

mlcp=N.MLCP(3,M,q)

def test_mlcp_enum_large():
    SO=N.SolverOptions(mlcp,N.SICONOS_MLCP_ENUM)
    N.mlcp_driver_init(mlcp, SO)
    info = N.mlcp_enum(mlcp, z, w, SO)
    N.mlcp_driver_reset(mlcp, SO)
    print "z = ", z
    print "w = ", w
    assert (linalg.norm(z-zsol) <= ztol)
    assert not info
mlcp =0
mlcp=N.MLCP()
N.mixedLinearComplementarity_newFromFilename(mlcp,"./data/diodeBridge_mlcp.dat")
#N.displayMLCP(mlcp)

def test_mlcp_enum_large_fromfile():
    SO=N.SolverOptions(mlcp,N.SICONOS_MLCP_ENUM)
    N.mlcp_driver_init(mlcp, SO)
    info = N.mlcp_enum(mlcp, z, w, SO)
    N.mlcp_driver_reset(mlcp, SO)
    print "z = ", z
    print "w = ", w
    assert (linalg.norm(z-zsol) <= ztol)
    assert not info
