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
import scipy.sparse as sp


def create_fcpb(dim, nc):
    mFC = np.zeros((dim, dim), dtype=np.float64, order="F")
    qFC = np.zeros(dim, dtype=np.float64)
    diag = np.ones(dim, dtype=np.float64)
    np.fill_diagonal(mFC, diag)

    qFC[:9] = [-1, 1, 3, -1, 1, 3, -1, 1, 3]

    problem = sn.FrictionContactProblem(3, nc, mFC, qFC)
    problem.M[1, 1] = 18
    mFC[3, 3] = 12
    assert mFC[1, 1] == 18
    assert problem.M[3, 3] == 12
    assert np.allclose(mFC, problem.M)
    assert np.allclose(qFC, problem.q)
    assert np.may_share_memory(mFC, problem.M)
    assert np.may_share_memory(qFC, problem.q)
    return problem


def create_fcpb2(dim, nc):
    qFC = np.ones(dim, dtype=np.float64)
    return sn.FrictionContactProblem(
        3, nc, np.ones((dim, dim), dtype=np.float64, order="F"), qFC
    )


def create_fcpb_sparse(dim, nc):

    data = np.array(
        [1.0, 0.5, -0.2, 1.0, 0.5, -0.2, 1.0, 0.5, -0.2, 1.0, 0.5, -0.2, 1.0, 0.5, 1.0],
        dtype=np.float64,
    )
    indices = np.array([0, 1, 2, 1, 2, 3, 2, 3, 4, 3, 4, 5, 4, 5, 5], dtype=np.int64)
    indptr = np.array([0, 3, 6, 9, 12, 14, 15], dtype=np.int64)
    mFC = sp.csc_matrix((data, indices, indptr), shape=(6, 6), dtype=np.float64)
    qFC = np.ones(dim, dtype=np.float64)

    problem = sn.FrictionContactProblem(3, nc, mFC, qFC)
    # Some prints to tests access and memory corruption
    print(problem.M)
    print("shape :", mFC.shape)
    print("indices :", mFC.indices)
    print("content :", mFC.data)
    print("indptr :", mFC.indptr)
    print(problem.M)
    problem.M[1, 1] = 18
    mFC[3, 3] = 12
    assert mFC[1, 1] == 18
    assert problem.M[3, 3] == 12
    assert np.shares_memory(mFC.data, problem.M.data)
    assert np.shares_memory(mFC.indices, problem.M.indices)
    assert np.shares_memory(mFC.indptr, problem.M.indptr)
    assert np.may_share_memory(qFC, problem.q)
    return problem


def test_friction_contact_build():

    nc = 3
    dim = 3 * nc
    problem = create_fcpb(dim, nc)
    problem.mu = [0.3, 0.2, 0.3]
    assert problem.q.size == dim
    assert problem.M.shape == (dim, dim)
    print(problem.M)  # just to check if it's still alive
    print(problem.q)  # just to check if it's still alive
    del problem

    pb2 = create_fcpb2(dim, nc)
    assert pb2.q.size == dim
    assert pb2.M.shape == (dim, dim)
    print(pb2.M)  # just to check if it's still alive
    print(pb2.q)  # just to check if it's still alive
    del pb2


def test_friction_contact_sparse_build():

    nc = 3
    dim = 3 * nc
    problem = create_fcpb_sparse(dim, nc)
    problem.mu = [0.3, 0.2, 0.3]
    assert problem.q.size == dim
    assert problem.M.shape == (dim, dim)
    print(problem.M)  # just to check if it's still alive
    print(problem.q)  # just to check if it's still alive
    del problem


def create_rfcpb(dim, nc):
    mFC = np.zeros((dim, dim), dtype=np.float64, order="F")
    qFC = np.zeros(dim, dtype=np.float64)
    diag = np.ones(10, dtype=np.float64)
    np.fill_diagonal(mFC, diag)

    qFC[...] = [-1, 1, 3, -1, 1, 3, -1, 1, 3, -1]

    problem = sn.RollingFrictionContactProblem(5, nc, mFC, qFC)

    problem.mu = [0.3, 0.2]
    assert problem.q.size == 10
    problem.M[1, 1] = 18
    mFC[3, 3] = 12
    assert problem.M.shape == (10, 10)
    assert mFC[1, 1] == 18
    assert problem.M[3, 3] == 12
    assert np.allclose(mFC, problem.M)
    assert np.allclose(qFC, problem.q)
    return problem


def create_rfcpb_sparse(dim, nc):

    data = np.array(
        [1.0, 0.5, -0.2, 1.0, 0.5, -0.2, 1.0, 0.5, -0.2, 1.0, 0.5, -0.2, 1.0, 0.5, 1.0],
        dtype=np.float64,
    )
    indices = np.array([0, 1, 2, 1, 2, 3, 2, 3, 4, 3, 4, 5, 4, 5, 5], dtype=np.int64)
    indptr = np.array([0, 3, 6, 9, 12, 14, 15], dtype=np.int64)
    mFC = sp.csc_matrix((data, indices, indptr), shape=(6, 6), dtype=np.float64)
    qFC = np.ones(dim, dtype=np.float64)

    problem = sn.RollingFrictionContactProblem(5, nc, mFC, qFC)
    # Some prints to tests access and memory corruption
    print(problem.M)
    print("shape :", mFC.shape)
    print("indices :", mFC.indices)
    print("content :", mFC.data)
    print("indptr :", mFC.indptr)
    print(problem.M)
    problem.M[1, 1] = 18
    mFC[3, 3] = 12
    assert mFC[1, 1] == 18
    assert problem.M[3, 3] == 12
    assert np.shares_memory(mFC.data, problem.M.data)
    assert np.shares_memory(mFC.indices, problem.M.indices)
    assert np.shares_memory(mFC.indptr, problem.M.indptr)
    assert np.may_share_memory(qFC, problem.q)
    return problem


def test_rolling_friction_contact_build():

    nc = 2
    dim = 5 * nc
    pb = create_rfcpb(dim, nc)
    print(pb.M)  # just to check if it's still alive
    print(pb.q)  # just to check if it's still alive
    del pb


def test_rolling_friction_contact_sparse_build():

    nc = 2
    dim = 5 * nc
    pb = create_rfcpb_sparse(dim, nc)
    print(pb.M)  # just to check if it's still alive
    print(pb.q)  # just to check if it's still alive
    del pb


def create_gfcpb(dim, nc, n):
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

    problem = sn.GlobalFrictionContactProblem(3, nc, mFC, H, qFC)

    problem.M[1, 1] = 4
    problem.H[1, 1] = 8
    H[0, 2] = 199
    assert np.allclose(mFC, problem.M)
    assert np.allclose(H, problem.H)
    assert np.allclose(qFC, problem.q)
    return problem


def create_gfcpb_sparse(dim, nc, n):
    nnz = 15
    data = np.random.random(nnz)  # dtype=np.float64)
    indices = np.array([0, 1, 2, 1, 2, 3, 2, 3, 4, 3, 4, 5, 4, 5, 5], dtype=np.int64)
    indptr = np.array([0, 3, 6, 9, 12, 14, 15], dtype=np.int64)
    mFC = sp.csc_matrix((data, indices, indptr), shape=(n, n), dtype=np.float64)
    qFC = np.ones(dim, dtype=np.float64)

    d2 = np.random.random(nnz)  # , dtype=np.float64)
    i2 = np.array([0, 1, 2, 1, 2, 3, 2, 3, 4, 3, 4, 5, 4, 5, 5], dtype=np.int64)
    iptr2 = np.array([0, 3, 6, 9, 12, 14, 15], dtype=np.int64)
    H = sp.csc_matrix((d2, i2, iptr2), shape=(n, dim * nc), dtype=np.float64)
    problem = sn.GlobalFrictionContactProblem(3, nc, mFC, H, qFC)
    # Some prints to tests access and memory corruption
    print(problem.M)
    print("shape :", mFC.shape)
    print("indices :", mFC.indices)
    print("content :", mFC.data)
    print("indptr :", mFC.indptr)
    print(problem.M)
    problem.M[1, 1] = 18
    mFC[3, 3] = 12
    assert mFC[1, 1] == 18
    assert problem.M[3, 3] == 12
    assert np.shares_memory(mFC.data, problem.M.data)
    assert np.shares_memory(mFC.indices, problem.M.indices)
    assert np.shares_memory(mFC.indptr, problem.M.indptr)
    assert np.may_share_memory(qFC, problem.q)
    return problem


def test_global_friction_contact_build():
    nc = 3
    dim = 3 * nc
    n = 7
    problem = create_gfcpb(dim, nc, n)
    print(problem.M)  # just to check if it's still alive
    print(problem.q)  # just to check if it's still alive
    del problem


def create_grfcpb(dim, nc, n):
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

    problem = sn.GlobalRollingFrictionContactProblem(5, nc, mFC, H, qFC)

    problem.M[1, 1] = 4
    problem.H[1, 1] = 8
    H[0, 2] = 199
    assert np.allclose(mFC, problem.M)
    assert np.allclose(H, problem.H)
    assert np.allclose(qFC, problem.q)
    return problem


def test_global_rolling_friction_contact_build():

    nc = 3
    dim = 5 * nc
    n = 7
    problem = create_grfcpb(dim, nc, n)
    print(problem.M)  # just to check if it's still alive
    print(problem.H)  # just to check if it's still alive
    print(problem.q)  # just to check if it's still alive
    del problem
