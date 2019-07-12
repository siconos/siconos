#!/usr/bin/env python
# Copyright (C) 2005, 2018 by INRIA

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
    mcp=SN.MixedComplementarityProblem(1, 1, mcp_function, mcp_Nablafunction)



def test_mcp_newton_FB_FBLSA():
    mcp = SN.MixedComplementarityProblem(0, 2, mcp_function, mcp_Nablafunction)
    z = np.array([0., 0.])
    w = np.array([0., 0.])

    SO = SN.SolverOptions(mcp, SN.SICONOS_MCP_NEWTON_FB_FBLSA)
    info = SN.mcp_newton_FB_FBLSA(mcp, z, w, SO)
    print("z = ", z)
    print("w = ", w)
    assert (np.linalg.norm(z-zsol) <= ztol)
    assert not info

def test_mcp_newton_min_FBLSA():
    mcp = SN.MixedComplementarityProblem(0, 2, mcp_function, mcp_Nablafunction)
    z = np.array([0., 0.])
    w = np.array([0., 0.])

    SO = SN.SolverOptions(mcp, SN.SICONOS_MCP_NEWTON_MIN_FBLSA)
    info = SN.mcp_newton_min_FBLSA(mcp, z, w, SO)
    print("z = ", z)
    print("w = ", w)
    assert (np.linalg.norm(z-zsol) <= ztol)
    assert not info

n=10
def build_problem(n):
    M = np.zeros((n,n))
    q = np.zeros(n)
    
    for i in range(n):
        q[i] = -i-5
        M[i,i] =2
        if i < n-1 :
            M[i,i+1] =1
        if i > 0 :
            M[i,i-1] =1
            
    return M,q



def mcp_function_2(n, z, F):
    M,q = build_problem(n)
    # F=np.dot(M,z) + q  pointer assignment is not working 
    F[:] = np.dot(M,z) + q
    return

def mcp_Nablafunction_2(n, z, nablaF):    
    M,q = build_problem(n)
    # nablaF= M pointer assignment is not working 
    nablaF[:] = M
    return 

def test_mcp_newton_FB_FBLSA_2():
    mcp = SN.MixedComplementarityProblem(n-3, 3, mcp_function_2, mcp_Nablafunction_2)
    z = np.zeros(n)
    w = np.zeros(n)
    SO = SN.SolverOptions(mcp, SN.SICONOS_MCP_NEWTON_FB_FBLSA)
    info = SN.mcp_newton_FB_FBLSA(mcp, z, w, SO)
    print("z = ", z)
    print("w = ", w)
    #assert (np.linalg.norm(z-zsol) <= ztol)
    assert not info

def test_mcp_newton_min_FBLSA_2():
    mcp = SN.MixedComplementarityProblem(n-3, 3, mcp_function_2, mcp_Nablafunction_2)
    z = np.zeros(n)
    w = np.zeros(n)
    SO = SN.SolverOptions(mcp, SN.SICONOS_MCP_NEWTON_MIN_FBLSA)
    info = SN.mcp_newton_min_FBLSA(mcp, z, w, SO)
    print("z = ", z)
    print("w = ", w)
    #assert (np.linalg.norm(z-zsol) <= ztol)
    assert not info

    
if __name__ == "__main__":
    SN.numerics_set_verbose(3)
    test_mcp_newton_FB_FBLSA()
    test_mcp_newton_min_FBLSA()
    test_mcp_newton_FB_FBLSA_2()
    test_mcp_newton_min_FBLSA_2()
