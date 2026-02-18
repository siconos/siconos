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
import siconos.numerics as sn


NC = 1

M = np.eye(3 * NC)
q = np.array([-1.0, 1.0, 3.0])
eta = np.array([0.1])
mu = np.array([0.1])
stress = np.array([0.0, 0.0, 0.0])
plastic_strain_rate = np.array([0.0, 0.0, 0.0])

sn.numerics_set_verbose(0)
MC = sn.MohrCoulomb2DProblem(3, M, q, eta, mu)
sn.mohrCoulomb2D_display(MC)

FC = sn.FrictionContactProblem(3, M, q, eta)


def solve(problem, solver, options):
    """Solve problem for a given solver"""
    stress[...] = 0.0
    plastic_strain_rate[...] = 0.0
    r = solver(problem, stress, plastic_strain_rate, options)
    size = stress.shape[0] // 3
    max_size = 3
    size_display = min(size, max_size)
    for i in range(size_display):
        for d in range(3):
            print(
                f"plastic_strain_rate[{3*i+d}] = {plastic_strain_rate[3*i+d]} \t\tstress[{3*i+d}] = {stress[3*i+d]}"
            )

    assert options.dparam[1] < options.dparam[0]
    assert not r


def test_mc2dnsgs():
    """Non-smooth Gauss Seidel, default"""
    SO = sn.SolverOptions(sn.MOHR_COULOMB_2D_NSGS)
    solve(MC, sn.mc2d_nsgs, SO)


def test_fc3dnsgs():
    """Non-smooth Gauss Seidel, default"""
    SO = sn.SolverOptions(sn.solver_ids.SICONOS_FRICTION_3D_NSGS)
    sn.solver_options_update_internal(
        SO, 0, sn.solver_ids.SICONOS_ONECONE_NSN
    )
    solve(FC, sn.fc3d_nsgs, SO)


def test_driver():
    """Non-smooth Gauss Seidel, default"""
    SO = sn.SolverOptions(sn.MOHR_COULOMB_2D_NSGS)
    r = sn.mc2d_driver(MC, stress, plastic_strain_rate, SO)


if __name__ == "__main__":
    test_driver()
    test_mc2dnsgs()
    test_fc3dnsgs()
