# Copyright 2022 INRIA
import numpy as np

import siconos.numerics as sn


def mcp_function(z):
    M = np.array([[2.0, 1.0], [1.0, 2.0]])

    q = np.array([-5.0, -6.0])
    return np.dot(M, z) + q


def mcp_Nablafunction(z):
    M = np.array([[2.0, 1.0], [1.0, 2.0]])
    return M


# solution
zsol = np.array([4.0 / 3.0, 7.0 / 3.0])
wsol = np.array([0.0, 0.0])

ztol = 1e-8


def test_new():
    return sn.MCP_old(1, 1, mcp_function, mcp_Nablafunction)


def test_mcp_FB():
    mcp = sn.MCP_old(1, 1, mcp_function, mcp_Nablafunction)
    z = np.array([0.0, 0.0])
    w = np.array([0.0, 0.0])

    options = sn.SolverOptions(sn.SICONOS_MCP_OLD_FB)
    sn.mcp_old_driver_init(mcp, options)
    info = sn.mcp_old_FischerBurmeister(mcp, z, w, options)
    sn.mcp_old_driver_reset(mcp, options)
    print("z = ", z)
    print("w = ", w)
    assert np.linalg.norm(z - zsol) <= ztol
    assert not info


n = 10


def build_problem(n):
    M = np.zeros((n, n), dtype=np.float64)
    q = np.zeros(n, dtype=np.float64)

    for i in range(n):
        q[i] = -i - 5
        M[i, i] = 2
        if i < n - 1:
            M[i, i + 1] = 1
        if i > 0:
            M[i, i - 1] = 1

    return M, q


def mcp_function_2(z):
    M, q = build_problem(n)
    return np.dot(M, z) + q


def mcp_Nablafunction_2(z):
    M, q = build_problem(n)
    return M


def test_mcp_FB_2():
    mcp = sn.MCP_old(n - 3, 3, mcp_function_2, mcp_Nablafunction_2)
    z = np.zeros(n)
    w = np.zeros(n)

    options = sn.SolverOptions(sn.SICONOS_MCP_OLD_FB)
    sn.mcp_old_driver_init(mcp, options)
    info = sn.mcp_old_FischerBurmeister(mcp, z, w, options)
    sn.mcp_old_driver_reset(mcp, options)
    print("z = ", z)
    print("w = ", w)
    assert not info


if __name__ == "__main__":
    sn.numerics_set_verbose(3)
    test_mcp_FB()
    test_mcp_FB_2()
