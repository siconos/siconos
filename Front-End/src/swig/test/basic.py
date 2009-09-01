#!/usr/bin/env python

# basic interface

from numpy import *

from Siconos import *


# Murty88, p2
M = array([[2., 1.],
           [1., 2.]])

Q = array([-5., -6.])

Z = array([0., 0.])

W = array([0., 0.])

# solution
Zsol = array([4./3., 7./3.])
Wsol = array([0. , 0.])

SO = Solver_Options("", [1000, 0], [0.001, 0.],1)

setNumericsVerbose(2)

#raw_input("Press Enter to terminate")

info  = lcp_pgs(LCP(M,Q),Z,W,SO)

print Z-Zsol
print W

info  = lcp_qp(LCP(M,Q),Z,W,SO)

print Z-Zsol
print W

info = lcp_lexicolemke(LCP(M,Q), Z, W, SO)

print Z-Zsol
print W

