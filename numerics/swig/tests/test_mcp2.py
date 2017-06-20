#!/usr/bin/env python
# Copyright (C) 2005, 2014 by INRIA

import numpy as np

# import siconos.numerics * fails with py.test!
import siconos.numerics as SN

def mcp_function(n, z, F):
    M = np.array([[2., 1.],
               [1., 2.]])

    q = np.array([-5., -6.])
    F[:] = np.dot(M,z) + q
    pass

def mcp_Nablafunction (n, z, nabla_F):
    M = np.array([[2., 1.],
               [1., 2.]])
    nabla_F[:] = M
    pass

# solution
zsol = np.array([4./3., 7./3.])
wsol = np.array([0. , 0.])

# problem
#mcp=N.MCP(1,1,mcp_function,mcp_Nablafunction)

ztol = 1e-8


def test_new():
    mcp=SN.MixedComplementarityProblem2(1, 1, mcp_function, mcp_Nablafunction)



def test_mcp_newton_FBLSA():
    mcp = SN.MixedComplementarityProblem2(0, 2, mcp_function, mcp_Nablafunction)
    z = np.array([0., 0.])
    w = np.array([0., 0.])

    SO = SN.SolverOptions(mcp, SN.SICONOS_MCP_NEWTON_FBLSA)
    info = SN.mcp_newton_FBLSA(mcp, z, w, SO)
    #print("z = ", z)
    #print("w = ", w)
    assert (np.linalg.norm(z-zsol) <= ztol)
    assert not info

def test_mcp_newton_minFBLSA():
    mcp = SN.MixedComplementarityProblem2(0, 2, mcp_function, mcp_Nablafunction)
    z = np.array([0., 0.])
    w = np.array([0., 0.])

    SO = SN.SolverOptions(mcp, SN.SICONOS_MCP_NEWTON_MINFBLSA)
    info = SN.mcp_newton_minFBLSA(mcp, z, w, SO)
    #print("z = ", z)
    #print("w = ", w)
    assert (np.linalg.norm(z-zsol) <= ztol)
    assert not info
