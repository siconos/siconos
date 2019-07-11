# Copyright (C) 2005, 2018 by INRIA
#!/usr/bin/env python
import numpy as np

import siconos.numerics as N

def mcp_function(z):
    M = np.array([[2., 1.],
                  [1., 2.]])

    q = np.array([-5., -6.])
    return np.dot(M,z) + q

def mcp_Nablafunction(z):
    M = np.array([[2., 1.],
                  [1., 2.]])
    return M

# solution
zsol = np.array([4./3., 7./3.])
wsol = np.array([0. , 0.])

# problem
#mcp=N.MCP(1,1,mcp_function,mcp_Nablafunction)

ztol = 1e-8


def test_new():
    mcp=N.MCP(1, 1, mcp_function, mcp_Nablafunction)

def test_mcp_FB():
    mcp=N.MCP(1,1,mcp_function,mcp_Nablafunction)
    z = np.array([0., 0.])
    w = np.array([0., 0.])

    SO=N.SolverOptions(mcp,N.SICONOS_MCP_FB)
    N.mcp_driver_init(mcp, SO)
    info = N.mcp_FischerBurmeister(mcp, z, w, SO)
    N.mcp_driver_reset(mcp, SO)
    print("z = ", z)
    print("w = ", w)
    assert (np.linalg.norm(z-zsol) <= ztol)
    assert not info

def test_mcp_newton_FBLSA():
    mcp=N.MCP(1,1,mcp_function,mcp_Nablafunction)
    z = np.array([0., 0.])
    w = np.array([0., 0.])

    SO=N.SolverOptions(mcp,N.SICONOS_MCP_FB)
    N.mcp_driver_init(mcp, SO)
    info = N.mcp_FischerBurmeister(mcp, z, w, SO)
    N.mcp_driver_reset(mcp, SO)
    print("z = ", z)
    print("w = ", w)
    assert (np.linalg.norm(z-zsol) <= ztol)
    assert not info


if __name__ == "__main__":
    N.numerics_set_verbose(3)
    test_mcp_FB()
    test_mcp_newton_FBLSA()
