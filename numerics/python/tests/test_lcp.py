# Siconos is a program dedicated to modeling, simulation and control
# of non smooth dynamical systems.
#
# Copyright 2026 INRIA.
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
# http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.
#

import siconos.numerics as sn
import numpy as np

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
    sn.solver_ids.SICONOS_LCP_PGS,
    sn.solver_ids.SICONOS_LCP_QP,
    sn.solver_ids.SICONOS_LCP_LEMKE,
    sn.solver_ids.SICONOS_LCP_ENUM,
]


def lcp_generic(id, z, w):

    options = sn.SolverOptions(id)
    info = sn.linearComplementarity_driver(lcp, z, w, options)
    print(" iter =", options.iparam[sn.params.SICONOS_IPARAM_ITER_DONE])
    print(" error=", options.dparam[sn.params.SICONOS_DPARAM_RESIDU])
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
