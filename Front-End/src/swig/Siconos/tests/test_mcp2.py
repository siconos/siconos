# Copyright (C) 2005, 2014 by INRIA
#!/usr/bin/env python

from numpy import *

# import Siconos.Numerics * fails with py.test!
import Siconos.Numerics as SN

def mcp_function(n1, n2, z, F) :
    M = array([[2., 1.],
               [1., 2.]])

    q = array([-5., -6.])
    F[:] = dot(M,z) + q
    pass

def mcp_Nablafunction (n1, n2, z, nabla_F) :
    M = array([[2., 1.],
               [1., 2.]])
    nabla_F[:] = M
    pass

# solution
zsol = array([4./3., 7./3.])
wsol = array([0. , 0.])

# problem
#mcp=N.MCP(1,1,mcp_function,mcp_Nablafunction)

ztol = 1e-8


def test_new():
    mcp=SN.MixedComplementarityProblem2(1, 1, mcp_function, mcp_Nablafunction)



def test_mcp_newton_FBLSA():
    mcp = SN.MixedComplementarityProblem2(0, 2, mcp_function, mcp_Nablafunction)
    z = array([0., 0.])
    w = array([0., 0.])

    SO = SN.SolverOptions(mcp, SN.SICONOS_MCP_NEWTON_FBLSA)
    info = SN.mcp_newton_FBLSA(mcp, z, w, SO)
    print("z = ", z)
    print("w = ", w)
    assert (linalg.norm(z-zsol) <= ztol)
    assert not info

def test_mcp_newton_minFBLSA():
    mcp = SN.MixedComplementarityProblem2(0, 2, mcp_function, mcp_Nablafunction)
    z = array([0., 0.])
    w = array([0., 0.])

    SO = SN.SolverOptions(mcp, SN.SICONOS_MCP_NEWTON_MINFBLSA)
    info = SN.mcp_newton_minFBLSA(mcp, z, w, SO)
    print("z = ", z)
    print("w = ", w)
    assert (linalg.norm(z-zsol) <= ztol)
    assert not info
