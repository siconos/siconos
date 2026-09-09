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
theta = np.array([0.1])
stress = np.array([0.0, 0.0, 0.0])
plastic_strain_rate = np.array([0.0, 0.0, 0.0])

sn.numerics_set_verbose(0)


def test_plasticity_problem_creation():
    """Test creating a PlasticityProblem"""
    P2D = sn.PlasticityProblem_new(3, M, q, eta, theta)
    assert P2D.dimension == 3
    assert P2D.numberOfCones == 1
    sn.plasticity_display(P2D)


def test_plasticity_problem_2d_alias():
    """Test that Plasticity2DProblem is an alias for PlasticityProblem"""
    P2D = sn.Plasticity2DProblem_new(3, M, q, eta, theta)
    assert P2D.dimension == 3
    assert P2D.numberOfCones == 1


def test_solver_ids():
    """Test that solver IDs are exported"""
    assert sn.PLASTICITY_2D_NSGS == 20000
    assert sn.MOHR_COULOMB_2D_NSGS == 20000  # backward compatibility
    assert sn.SICONOS_MC2D_NSGS == 20000  # backward compatibility
    assert sn.PLASTICITY_2D_ONECONE_NSN == 20050


def test_solver_options_create():
    """Test creating solver options"""
    options = sn.SolverOptions(sn.solver_ids.PLASTICITY_2D_NSGS)
    assert options is not None


if __name__ == "__main__":
    test_plasticity_problem_creation()
    test_plasticity_problem_2d_alias()
    test_solver_ids()
    test_solver_options_create()
    print("All tests passed!")
