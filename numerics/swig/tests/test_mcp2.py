#!/usr/bin/env python
# Copyright 2022 INRIA

import numpy as np

# import siconos.numerics * fails with py.test!
import siconos.numerics as sn

# solution
zsol = np.array([4.0 / 3.0, 7.0 / 3.0])
wsol = np.array([0.0, 0.0])
ztol = 1e-8


def mcp_function(n, z, F):
    M = np.array([[2.0, 1.0], [1.0, 2.0]])

    q = np.array([-5.0, -6.0])
    F[:] = np.dot(M, z) + q
    pass


def mcp_Nablafunction(n, z, nabla_F):
    M = np.array([[2.0, 1.0], [1.0, 2.0]])
    nabla_F[:] = M
    pass


def test_new():
    mcp = sn.MCP(1, 1, mcp_function, mcp_Nablafunction)
    z = np.array([0.0, 0.0])
    w = np.array([0.0, 0.0])

    SO = sn.SolverOptions(sn.SICONOS_MCP_NEWTON_FB_FBLSA)
    info = sn.mcp_newton_FB_FBLSA(mcp, z, w, SO)
    print("z = ", z)
    print("w = ", w)
    assert np.linalg.norm(z - zsol) <= ztol
    assert not info


def test_mcp_newton_FB_FBLSA():
    mcp = sn.MCP(0, 2, mcp_function, mcp_Nablafunction)
    z = np.array([0.0, 0.0])
    w = np.array([0.0, 0.0])

    SO = sn.SolverOptions(sn.SICONOS_MCP_NEWTON_FB_FBLSA)
    info = sn.mcp_newton_FB_FBLSA(mcp, z, w, SO)
    print("z = ", z)
    print("w = ", w)
    assert np.linalg.norm(z - zsol) <= ztol
    assert not info


def test_mcp_newton_min_FBLSA():
    mcp = sn.MCP(0, 2, mcp_function, mcp_Nablafunction)
    z = np.array([0.0, 0.0])
    w = np.array([0.0, 0.0])

    SO = sn.SolverOptions(sn.SICONOS_MCP_NEWTON_MIN_FBLSA)
    info = sn.mcp_newton_min_FBLSA(mcp, z, w, SO)
    print("z = ", z)
    print("w = ", w)
    assert np.linalg.norm(z - zsol) <= ztol
    assert not info


def build_problem(n):
    M = np.zeros((n, n))
    q = np.zeros(n)

    for i in range(n):
        q[i] = -i + 7
        M[i, i] = 2
        if i < n - 1:
            M[i, i + 1] = 1
        if i > 0:
            M[i, i - 1] = 1

    return M, q


def mcp_function_2(n, z, F):
    M, q = build_problem(n)
    # F=np.dot(M,z) + q  pointer assignment is not working
    F[:] = np.dot(M, z) + q
    return


def mcp_Nablafunction_2(n, z, nablaF):
    M, q = build_problem(n)
    # nablaF= M pointer assignment is not working
    nablaF[:] = M
    return


def test_mcp_newton_FB_FBLSA_2():
    n = 10
    mcp = sn.MCP(n - 5, 5, mcp_function_2, mcp_Nablafunction_2)
    z = np.zeros(n)
    w = np.zeros(n)
    SO = sn.SolverOptions(sn.SICONOS_MCP_NEWTON_FB_FBLSA)
    info = sn.mcp_newton_FB_FBLSA(mcp, z, w, SO)
    print("z = ", z)
    print("w = ", w)
    # assert (np.linalg.norm(z-zsol) <= ztol)
    assert not info


def test_mcp_newton_min_FBLSA_2():
    n = 10
    mcp = sn.MCP(n - 5, 5, mcp_function_2, mcp_Nablafunction_2)
    z = np.zeros(n)
    w = np.zeros(n)
    options = sn.SolverOptions(sn.SICONOS_MCP_NEWTON_MIN_FBLSA)
    options.iparam[
        sn.SICONOS_IPARAM_STOPPING_CRITERION
    ] = sn.SICONOS_STOPPING_CRITERION_RESIDU

    sn.solver_options_print(options)

    info = sn.mcp_newton_min_FBLSA(mcp, z, w, options)
    print("z = ", z)
    print("w = ", w)
    # assert (np.linalg.norm(z-zsol) <= ztol)
    assert not info


if __name__ == "__main__":
    sn.numerics_set_verbose(3)
    test_mcp_newton_FB_FBLSA()
    test_mcp_newton_min_FBLSA()
    test_mcp_newton_FB_FBLSA_2()
    test_mcp_newton_min_FBLSA_2()
