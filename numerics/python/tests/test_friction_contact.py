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

import numpy as np
import scipy.sparse as sp
import siconos.numerics as sn


def create_fcpb(dim, nc):
    mFC = np.zeros((dim, dim), dtype=np.float64, order="F")
    qFC = np.zeros(dim, dtype=np.float64)
    diag = np.ones(dim, dtype=np.float64)
    np.fill_diagonal(mFC, diag)

    qFC[:9] = [-1, 1, 3, -1, 1, 3, -1, 1, 3]
    mu = np.ones(nc, dtype=np.float64)
    problem = sn.FrictionContactProblem(3, nc, mFC, qFC, mu)
    problem.M()[1, 1] = 18
    mFC[3, 3] = 12
    assert mFC[1, 1] == 18
    assert problem.M()[3, 3] == 12
    assert np.allclose(mFC, problem.M())
    assert np.allclose(qFC, problem.q())
    assert np.may_share_memory(mFC, problem.M())
    assert np.may_share_memory(qFC, problem.q())
    assert np.allclose(problem.mu(), mu)
    return problem


def create_fcpb2(dim, nc):
    qFC = np.ones(dim, dtype=np.float64)
    mu = np.ones(nc, dtype=np.float64)
    return sn.FrictionContactProblem(
        3, nc, np.ones((dim, dim), dtype=np.float64, order="F"), qFC, mu
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
    mu = np.ones(nc, dtype=np.float64)
    problem = sn.FrictionContactProblem(3, nc, mFC, qFC, mu)
    # Some prints to tests access and memory corruption
    print(problem.M())
    print("shape :", mFC.shape)
    print("indices :", mFC.indices)
    print("content :", mFC.data)
    print("indptr :", mFC.indptr)
    print(problem.M())
    problem.M()[1, 1] = 18
    mFC[3, 3] = 12
    assert mFC[1, 1] == 18
    assert problem.M()[3, 3] == 12
    assert np.shares_memory(mFC.data, problem.M().data)
    assert np.may_share_memory(qFC, problem.q())
    assert np.allclose(problem.mu(), mu)
    return problem


def test_friction_contact_build():

    nc = 3
    dim = 3 * nc
    problem = create_fcpb(dim, nc)
    #    problem.mu() = [0.3, 0.2, 0.3]
    assert problem.q().size == dim
    assert problem.M().shape == (dim, dim)
    print(problem.M())  # just to check if it's still alive
    print(problem.q())  # just to check if it's still alive
    del problem

    pb2 = create_fcpb2(dim, nc)
    assert pb2.q().size == dim
    assert pb2.M().shape == (dim, dim)
    print(pb2.M())  # just to check if it's still alive
    print(pb2.q())  # just to check if it's still alive
    del pb2


def test_friction_contact_sparse_build():

    nc = 2
    dim = 3 * nc
    problem = create_fcpb_sparse(dim, nc)
    # problem.mu() = [0.3, 0.2, 0.3]
    assert problem.q().size == dim
    assert problem.M().shape == (dim, dim)
    print(problem.M())  # just to check if it's still alive
    print(problem.q())  # just to check if it's still alive
    del problem


def test_compute_rho():
    nc = 2
    dim = 3 * nc

    # Dense problem with non-trivial diagonal blocks
    M = np.eye(dim, dtype=np.float64, order="F")
    M[1, 1] = 4.0
    M[2, 2] = 9.0
    M[4, 4] = 2.0
    M[5, 5] = 5.0
    q = np.zeros(dim, dtype=np.float64)
    mu = np.ones(nc, dtype=np.float64)
    problem = sn.FrictionContactProblem(3, nc, M, q, mu)

    rho_strat = sn.params.SICONOS_FRICTION_3D_NSN_RHO_STRATEGY_enum

    # Constant strategy
    rho_const = problem.compute_rho(
        rho_strat.SICONOS_FRICTION_3D_NSN_FORMULATION_RHO_STRATEGY_CONSTANT, 0.5
    )
    assert rho_const.shape == (nc, 3)
    assert np.allclose(rho_const, 0.5)

    # Split spectral norm strategy
    rho_split = problem.compute_rho(
        rho_strat.SICONOS_FRICTION_3D_NSN_FORMULATION_RHO_STRATEGY_SPLIT_SPECTRAL_NORM
    )
    assert rho_split.shape == (nc, 3)
    assert np.allclose(rho_split[:, 0], 1.0)
    assert rho_split[0, 1] > 0.0
    assert rho_split[0, 2] > 0.0

    # Split spectral norm cond strategy
    rho_split_cond = problem.compute_rho(
        rho_strat.SICONOS_FRICTION_3D_NSN_FORMULATION_RHO_STRATEGY_SPLIT_SPECTRAL_NORM_COND
    )
    assert rho_split_cond.shape == (nc, 3)
    assert np.allclose(rho_split_cond[:, 0], 1.0)
    assert rho_split_cond[0, 1] > 0.0
    assert rho_split_cond[0, 2] > 0.0

    # Spectral norm strategy
    rho_spec = problem.compute_rho(
        rho_strat.SICONOS_FRICTION_3D_NSN_FORMULATION_RHO_STRATEGY_SPECTRAL_NORM
    )
    assert rho_spec.shape == (nc, 3)
    assert np.allclose(rho_spec[:, 0], 1.0)

    # Sparse problem
    data = np.array(
        [
            1.0,
            1.0,
            1.0,
            1.0,
            4.0,
            1.0,
            1.0,
            9.0,
            1.0,
            1.0,
            1.0,
            1.0,
            1.0,
            2.0,
            1.0,
            1.0,
            5.0,
            1.0,
        ],
        dtype=np.float64,
    )
    indices = np.array(
        [0, 1, 2, 0, 1, 2, 0, 1, 2, 3, 4, 5, 3, 4, 5, 3, 4, 5],
        dtype=np.int64,
    )
    indptr = np.array([0, 3, 6, 9, 12, 15, 18], dtype=np.int64)
    M_sparse = sp.csc_matrix(
        (data, indices, indptr), shape=(dim, dim), dtype=np.float64
    )
    problem_sparse = sn.FrictionContactProblem(3, nc, M_sparse, q, mu)
    rho_sparse = problem_sparse.compute_rho(
        rho_strat.SICONOS_FRICTION_3D_NSN_FORMULATION_RHO_STRATEGY_SPLIT_SPECTRAL_NORM
    )
    assert rho_sparse.shape == (nc, 3)
    assert np.allclose(rho_sparse[:, 0], 1.0)
    assert rho_sparse[0, 1] > 0.0
    assert rho_sparse[0, 2] > 0.0


def create_rfcpb(dim, nc):
    mFC = np.zeros((dim, dim), dtype=np.float64, order="F")
    qFC = np.zeros(dim, dtype=np.float64)
    diag = np.ones(10, dtype=np.float64)
    np.fill_diagonal(mFC, diag)

    qFC[...] = [-1, 1, 3, -1, 1, 3, -1, 1, 3, -1]
    mu = np.ones(nc, dtype=np.float64)
    mu_r = 2.0 * np.ones(nc, dtype=np.float64)
    problem = sn.RollingFrictionContactProblem(5, nc, mFC, qFC, mu, mu_r)

    # problem.mu() = [0.3, 0.2]
    assert problem.q().size == 10
    problem.M()[1, 1] = 18
    mFC[3, 3] = 12
    assert problem.M().shape == (10, 10)
    assert mFC[1, 1] == 18
    assert problem.M()[3, 3] == 12
    assert np.allclose(mFC, problem.M())
    assert np.allclose(qFC, problem.q())
    assert np.allclose(problem.mu(), mu)
    assert np.allclose(problem.mu_r(), mu_r)
    return problem


def create_rfcpb_sparse(dim, nc):

    nnz = 15
    mFC = sp.random(dim, dim, density=nnz / (dim * dim), format="csc", dtype=np.float64)
    qFC = np.ones(dim, dtype=np.float64)
    mu = np.ones(nc, dtype=np.float64)
    mu_r = 2 * np.ones(nc, dtype=np.float64)

    problem = sn.RollingFrictionContactProblem(5, nc, mFC, qFC, mu, mu_r)
    # Some prints to tests access and memory corruption
    print(problem.M())
    print("shape :", mFC.shape)
    print("indices :", mFC.indices)
    print("content :", mFC.data)
    print("indptr :", mFC.indptr)
    print(problem.M())
    problem.M().data[1] = 18
    mFC.data[2] = 12
    assert np.allclose(problem.M().data, mFC.data)
    assert np.shares_memory(mFC.data, problem.M().data)
    assert np.may_share_memory(qFC, problem.q())
    assert np.allclose(problem.mu(), mu)
    assert np.allclose(problem.mu_r(), mu_r)
    return problem


def test_rolling_friction_contact_build():

    nc = 2
    dim = 5 * nc
    pb = create_rfcpb(dim, nc)
    print(pb.M())  # just to check if it's still alive
    print(pb.q())  # just to check if it's still alive
    del pb


def test_rolling_friction_contact_sparse_build():

    nc = 2
    dim = 5 * nc
    pb = create_rfcpb_sparse(dim, nc)
    print(pb.M())  # just to check if it's still alive
    print(pb.q())  # just to check if it's still alive
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

    mu = np.ones(nc, dtype=np.float64)

    problem = sn.GlobalFrictionContactProblem(3, nc, mFC, H, qFC, mu)

    problem.M()[1, 1] = 4
    problem.H()[1, 1] = 8
    H[0, 2] = 199
    assert np.allclose(mFC, problem.M())
    assert np.allclose(H, problem.H())
    assert np.allclose(qFC, problem.q())
    assert np.allclose(problem.mu(), mu)
    return problem


def create_gfcpb_sparse(dim, nc, n):
    nnz = 15
    mFC = sp.random(n, n, density=nnz / (dim * dim), format="csc", dtype=np.float64)
    qFC = np.ones(n, dtype=np.float64)

    H = sp.random(n, dim, density=nnz / (dim * dim), format="csc", dtype=np.float64)
    mu = np.ones(nc, dtype=np.float64)
    problem = sn.GlobalFrictionContactProblem(3, nc, mFC, H, qFC, mu)
    # Some prints to tests access and memory corruption
    print("shape :", mFC.shape)
    print("indices :", mFC.indices)
    print("content :", mFC.data)
    print("indptr :", mFC.indptr)
    print(problem.M())
    problem.M().data[1] = 18
    mFC.data[2] = 12
    assert np.allclose(problem.M().data, mFC.data)
    assert np.shares_memory(mFC.data, problem.M().data)
    assert np.may_share_memory(qFC, problem.q())
    assert np.allclose(problem.mu(), mu)

    return problem


def test_global_friction_contact_build():
    nc = 3
    dim = 3 * nc
    n = 7
    problem = create_gfcpb(dim, nc, n)
    print(problem.M())  # just to check if it's still alive
    print(problem.q())  # just to check if it's still alive
    del problem


def test_global_friction_contact_sparse_build():
    nc = 3
    dim = 3 * nc
    n = 7
    problem = create_gfcpb_sparse(dim, nc, n)
    print(problem.M())  # just to check if it's still alive
    print(problem.q())  # just to check if it's still alive
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
    mu = np.ones(nc, dtype=np.float64)
    mu_r = 2 * np.ones(nc, dtype=np.float64)

    problem = sn.GlobalRollingFrictionContactProblem(5, nc, mFC, H, qFC, mu, mu_r)

    problem.M()[1, 1] = 4
    problem.H()[1, 1] = 8
    H[0, 2] = 199
    assert np.allclose(mFC, problem.M())
    assert np.allclose(H, problem.H())
    assert np.allclose(qFC, problem.q())
    assert np.allclose(problem.mu(), mu)
    assert np.allclose(problem.mu_r(), mu_r)

    return problem


def create_grfcpb_sparse(dim, nc, n):
    nnz = 15
    mFC = sp.random(n, n, density=nnz / (dim * dim), format="csc", dtype=np.float64)
    qFC = np.ones(n, dtype=np.float64)

    H = sp.random(n, dim, density=nnz / (dim * dim), format="csc", dtype=np.float64)
    mu = np.ones(nc, dtype=np.float64)
    mur = 2.0 * np.ones(nc, dtype=np.float64)
    problem = sn.GlobalRollingFrictionContactProblem(5, nc, mFC, H, qFC, mu, mur)
    # Some prints to tests access and memory corruption
    print("shape :", mFC.shape)
    print("indices :", mFC.indices)
    print("content :", mFC.data)
    print("indptr :", mFC.indptr)
    print(problem.M())
    problem.M().data[1] = 18
    mFC.data[2] = 12
    assert np.allclose(problem.M().data, mFC.data)
    assert np.shares_memory(mFC.data, problem.M().data)
    assert np.may_share_memory(qFC, problem.q())
    assert np.allclose(problem.mu(), mu)
    assert np.allclose(problem.mu_r(), mur)
    return problem


def test_global_rolling_friction_contact_build():

    nc = 2
    dim = 5 * nc
    n = 7
    problem = create_grfcpb(dim, nc, n)
    print(problem.M())  # just to check if it's still alive
    print(problem.H())  # just to check if it's still alive
    print(problem.q())  # just to check if it's still alive
    del problem


def test_global_rolling_friction_contact_sparse_build():
    nc = 2
    dim = 5 * nc
    n = 7
    problem = create_grfcpb_sparse(dim, nc, n)
    print(problem.M())  # just to check if it's still alive
    print(problem.q())  # just to check if it's still alive
    del problem
