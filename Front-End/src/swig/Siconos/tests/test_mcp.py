# Copyright (C) 2005, 2012 by INRIA
#!/usr/bin/env python

from numpy import *

# import Siconos.Numerics * fails with py.test!
import Siconos.Numerics as N

def mcp_function (z) :
    M = array([[2., 1.],
               [1., 2.]])
    
    q = array([-5., -6.])
    return dot(M,z) + q

def mcp_Nablafunction (z) : 
    print "in mcp_function"
    M = array([[2., 1.],
               [1., 2.]])
    return M

# solution
zsol = array([4./3., 7./3.])
wsol = array([0. , 0.])

# problem
#mcp=N.MCP(1,1,mcp_function,mcp_Nablafunction)

ztol = 1e-8


def test_new():
    mcp=N.MCP(1,1,mcp_function,mcp_Nablafunction)



def test_mcp_FB():
    pass
    # SO=N.SolverOptions(mlcp,N.SICONOS_MLCP_ENUM)
    # N.mlcp_driver_init(mlcp, SO)
    # info = N.mlcp_enum(mlcp, z, w, SO)
    # N.mlcp_driver_reset(mlcp, SO)
    # print "z = ", z
    # print "w = ", w
    # assert (linalg.norm(z-zsol) <= ztol)
    # assert not info
