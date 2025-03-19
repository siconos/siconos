# Siconos is a program dedicated to modeling, simulation and control
# of non smooth dynamical systems.
#
# Copyright 2025 INRIA.
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


def test_friction_contact_build():

    nc = 3
    dim = 3 * nc
    mFC = np.zeros((dim, dim), dtype=np.float64, order="F")
    qFC = np.zeros(dim, dtype=np.float64)
    diag = np.ones(9, dtype=np.float64)
    np.fill_diagonal(mFC, diag)

    qFC[...] = [-1, 1, 3, -1, 1, 3, -1, 1, 3]

    problem = sn.FrictionContactProblem(3, nc)

    problem.q[...] = qFC
    problem.M = mFC
    problem.mu = [0.3, 0.2, 0.3]
    assert problem.q.size == 9
    problem.M[1, 1] = 18
    mFC[3, 3] = 12
    assert problem.M.shape == (9, 9)
    assert mFC[1, 1] == 18
    assert problem.M[3, 3] == 12
    assert np.allclose(mFC, problem.M)
    assert np.allclose(qFC, problem.q)


def test_rolling_friction_contact_build():

    nc = 2
    dim = 5 * nc
    mFC = np.zeros((dim, dim), dtype=np.float64, order="F")
    qFC = np.zeros(dim, dtype=np.float64)
    diag = np.ones(10, dtype=np.float64)
    np.fill_diagonal(mFC, diag)

    qFC[...] = [-1, 1, 3, -1, 1, 3, -1, 1, 3, -1]

    problem = sn.RollingFrictionContactProblem(5, nc)

    problem.q[...] = qFC
    problem.M = mFC
    problem.mu = [0.3, 0.2]
    assert problem.q.size == 10
    problem.M[1, 1] = 18
    mFC[3, 3] = 12
    assert problem.M.shape == (10, 10)
    assert mFC[1, 1] == 18
    assert problem.M[3, 3] == 12
    assert np.allclose(mFC, problem.M)
    assert np.allclose(qFC, problem.q)


def test_global_friction_contact_build():
    nc = 3
    dim = 3 * nc
    n = 7
    mFC = np.zeros((n, n), dtype=np.float64, order="F")
    qFC = np.zeros(n, dtype=np.float64)
    diag = np.ones(n, dtype=np.float64)
    np.fill_diagonal(mFC, diag)
    mFC[0, 4] = 17

    H = np.zeros((n, dim), dtype=np.float64, order="F")
    H[0, 0] = 1
    H[1, 2] = 6
    H[3, 2] = 8

    qFC[...] = [-1, 1, 3, -1, 1, 3, -1]

    problem = sn.GlobalFrictionContactProblem(3, nc)

    problem.q = qFC
    problem.M = mFC
    problem.H = H
    problem.M[1, 1] = 4
    problem.H[1, 1] = 8
    H[0, 2] = 199
    assert np.allclose(mFC, problem.M)
    assert np.allclose(H, problem.H)
    assert np.allclose(qFC, problem.q)


def test_global_rolling_friction_contact_build():

    nc = 3
    dim = 5 * nc
    n = 7
    mFC = np.zeros((n, n), dtype=np.float64, order="F")
    qFC = np.zeros(n, dtype=np.float64)
    diag = np.ones(n, dtype=np.float64)
    np.fill_diagonal(mFC, diag)
    mFC[0, 4] = 17

    H = np.zeros((n, dim), dtype=np.float64, order="F")
    H[0, 0] = 1
    H[1, 2] = 6
    H[3, 2] = 8

    qFC[...] = [-1, 1, 3, -1, 1, 3, -1]

    problem = sn.GlobalRollingFrictionContactProblem(5, nc)

    problem.q = qFC
    problem.M = mFC
    problem.H = H
    problem.M[1, 1] = 4
    problem.H[1, 1] = 8
    H[0, 2] = 199
    assert np.allclose(mFC, problem.M)
    assert np.allclose(H, problem.H)
    assert np.allclose(qFC, problem.q)
