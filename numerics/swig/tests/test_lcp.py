#!/usr/bin/env python
import numpy as np

# import Siconos.Numerics * fails with py.test!
import siconos.numerics as N
import siconos
# basic interface
# Murty88, p2
M = np.array([[2., 1.],
           [1., 2.]])

q = np.array([-5., -6.])

z = np.array([0., 0.])

w = np.array([0., 0.])

# solution
zsol = np.array([4./3., 7./3.])
wsol = np.array([0., 0.])

# problem
lcp = N.LCP(M, q)

ztol = 1e-4

def test_lcp_pgs():
    SO=N.SolverOptions(lcp,N.SICONOS_LCP_PGS)
    info  = N.lcp_pgs(lcp,z,w,SO)
    print('pgs iter =', SO.iparam[1])
    print('pgs error=', SO.dparam[1])
    assert (np.linalg.norm(z-zsol) <= ztol)
    assert not info

if siconos.WITH_FORTRAN and siconos.WITH_QL0001:
    def test_lcp_qp():
        SO=N.SolverOptions(lcp,N.SICONOS_LCP_QP)
        info  = N.lcp_qp(lcp,z,w,SO)
        assert (np.linalg.norm(z-zsol) <= ztol)
        assert not info

def test_lcp_lexicolemke():
    SO=N.SolverOptions(lcp, N.SICONOS_LCP_LEMKE)
    info = N.lcp_lexicolemke(lcp, z, w, SO)
    print('lexicolemke iter =', SO.iparam[1])
    assert (np.linalg.norm(z-zsol) <= ztol)
    assert not info

def test_lcp_enum():
    SO=N.SolverOptions(lcp,N.SICONOS_LCP_ENUM)
    info = N.lcp_enum(lcp, z, w, SO)
    assert (np.linalg.norm(z-zsol) <= ztol)
    assert not info
