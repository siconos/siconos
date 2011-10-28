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
lcp=N.LCP(M,q)

ztol = 1e-4

def test_lcp_pgs():
    SO=N.SolverOptions(lcp,N.SICONOS_LCP_PGS)
    info  = N.lcp_pgs(lcp,z,w,SO)
    print 'pgs iter =', SO.iparam[1]
    print 'pgs error=', SO.dparam[1]
    assert (linalg.norm(z-zsol) <= ztol)
    assert not info

def test_lcp_qp():
    SO=N.SolverOptions(lcp,N.SICONOS_LCP_QP)
    info  = N.lcp_qp(lcp,z,w,SO)
    assert (linalg.norm(z-zsol) <= ztol)
    assert not info

def test_lcp_lexicolemke():
    SO=N.SolverOptions(lcp,N.SICONOS_LCP_LEMKE)
    info = N.lcp_lexicolemke(lcp, z, w, SO)
    print 'lexicolemke iter =', SO.iparam[1]
    assert (linalg.norm(z-zsol) <= ztol)
    assert not info

def test_lcp_enum():
    SO=N.SolverOptions(lcp,N.SICONOS_LCP_ENUM)
    info = N.lcp_enum(lcp, z, w, SO)
    assert (linalg.norm(z-zsol) <= ztol)
    assert not info

