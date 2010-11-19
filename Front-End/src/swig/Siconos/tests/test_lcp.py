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


def test_lcp_pgs():
    SO=N.SolverOptions()
    N.linearComplementarity_setDefaultSolverOptions(lcp,SO,N.SICONOS_LCP_PGS)
    info  = N.lcp_pgs(lcp,z,w,SO)
    print z-zsol
    print w
    assert (linalg.norm(z-zsol) <= 0.01)
    assert not info

def test_lcp_qp():
    SO=N.SolverOptions()
    N.linearComplementarity_setDefaultSolverOptions(lcp,SO,N.SICONOS_LCP_QP)
    info  = N.lcp_qp(lcp,z,w,SO)

    print z-zsol
    print w
    assert (linalg.norm(z-zsol) <= 0.01)
    assert not info

def test_lcp_lexicolemke():
    SO=N.SolverOptions()
    N.linearComplementarity_setDefaultSolverOptions(lcp,SO,N.SICONOS_LCP_LEMKE)
    info = N.lcp_lexicolemke(lcp, z, w, SO)

    print z-zsol
    print w
    assert (linalg.norm(z-zsol) <= 0.01)
    assert not info


