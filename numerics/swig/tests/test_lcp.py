#!/usr/bin/env python
import numpy as np

# import Siconos.Numerics * fails with py.test!
import siconos.numerics as sn
import siconos

# basic interface
# Murty88, p2
M = np.array([[2.0, 1.0], [1.0, 2.0]])

q = np.array([-5.0, -6.0])

# solution
zsol = np.array([4.0 / 3.0, 7.0 / 3.0])
wsol = np.array([0.0, 0.0])

# problem
lcp = sn.LCP(M, q)

ztol = 1e-4


solvers = [
    sn.SICONOS_LCP_PGS,
    sn.SICONOS_LCP_QP,
    sn.SICONOS_LCP_LEMKE,
    sn.SICONOS_LCP_ENUM,
]


def lcp_generic(id, z, w):

    options = sn.SolverOptions(id)
    info = sn.linearComplementarity_driver(lcp, z, w, options)
    print(" iter =", options.iparam[sn.SICONOS_IPARAM_ITER_DONE])
    print(" error=", options.dparam[sn.SICONOS_DPARAM_RESIDU])
    return info


def test_lcps():

    z = np.zeros((2,), np.float64)
    w = np.zeros_like(z)

    for id in solvers:
        info = lcp_generic(id, z, w)
        assert np.linalg.norm(z - zsol) <= ztol
        assert not info
        z[...] = w[...] = 0.0


if __name__ == "__main__":

    test_lcps()
