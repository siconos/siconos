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

import numpy as np
import siconos

# import siconos.numerics * fails with py.test!
import siconos.numerics as SN


def ncp_function(n, z, F):
    M = np.array([[2.0, 1.0], [1.0, 2.0]])

    q = np.array([-5.0, -6.0])
    F[:] = np.dot(M, z) + q


def ncp_Nablafunction(n, z, nabla_F):
    M = np.array([[2.0, 1.0], [1.0, 2.0]])
    nabla_F[:] = M


# solution
zsol = np.array([4.0 / 3.0, 7.0 / 3.0])
wsol = np.array([0.0, 0.0])

# problem
# ncp=N.NCP(1,1,ncp_function,ncp_Nablafunction)

ztol = 1e-8


def test_new():
    return SN.NCP(1, ncp_function, ncp_Nablafunction)


def test_ncp_newton_FBLSA():
    ncp = SN.NCP(2, ncp_function, ncp_Nablafunction)
    z = np.array([0.0, 0.0])
    w = np.array([0.0, 0.0])

    SO = SN.SolverOptions(SN.SICONOS_NCP_NEWTON_FB_FBLSA)
    info = SN.ncp_driver(ncp, z, w, SO)
    assert np.linalg.norm(z - zsol) <= ztol
    assert not info


def test_ncp_newton_minFBLSA():
    ncp = SN.NCP(2, ncp_function, ncp_Nablafunction)
    z = np.array([0.0, 0.0])
    w = np.array([0.0, 0.0])

    SO = SN.SolverOptions(SN.SICONOS_NCP_NEWTON_MIN_FBLSA)
    info = SN.ncp_driver(ncp, z, w, SO)
    # print("z = ", z)
    # print("w = ", w)
    assert np.linalg.norm(z - zsol) <= ztol
    assert not info


def test_ncp_path():
    ncp = SN.NCP(2, ncp_function, ncp_Nablafunction)
    z = np.array([0.0, 0.0])
    w = np.array([0.0, 0.0])

    SO = SN.SolverOptions(SN.SICONOS_NCP_PATH)
    info = SN.ncp_driver(ncp, z, w, SO)

    if siconos.WITH_PATHFERRIS:
        assert np.linalg.norm(z - zsol) <= ztol
        assert not info
        return
    else:
        assert info != 0

        try:
            SN.ncp_path(ncp, z, w, SO)
        except RuntimeError:
            pass


if __name__ == "__main__":
    SN.numerics_set_verbose(3)
    test_ncp_newton_FBLSA()
    test_new()
    test_ncp_newton_minFBLSA()
    test_ncp_path()
