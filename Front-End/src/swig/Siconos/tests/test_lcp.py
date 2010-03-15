#!/usr/bin/env python

from numpy import *
from Siconos.Numerics import *

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
lcp=LCP(M,q)


def test_lcp_pgs():
    SO=SolverOptions()
    linearComplementarity_setDefaultSolverOptions(lcp,SO,"PGS")
    info  = lcp_pgs(lcp,z,w,SO)
    print z-zsol
    print w
    assert (linalg.norm(z-zsol) <= 0.01)

def test_lcp_qp():
    SO=SolverOptions()
    linearComplementarity_setDefaultSolverOptions(lcp,SO,"QP")
    info  = lcp_qp(LCP(M,q),z,w,SO)

    print z-zsol
    print w
    assert (linalg.norm(z-zsol) <= 0.01)

def test_lcp_lexicolemke():
    SO=SolverOptions()
    linearComplementarity_setDefaultSolverOptions(lcp,SO,"Lemke")
    info = lcp_lexicolemke(LCP(M,q), z, w, SO)

    print z-zsol
    print w
    assert (linalg.norm(z-zsol) <= 0.01)


