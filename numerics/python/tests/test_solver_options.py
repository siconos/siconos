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
from siconos.numerics import solver_ids
from siconos.numerics import params as pnames


def test_solver_options_create():
    sid = solver_ids.SICONOS_FRICTION_3D_NSGS
    #sn.numerics_set_verbose(2)
    so = sn.SolverOptions(sid)
    so.iparam[pnames.SICONOS_IPARAM_MAX_ITER] = 1000
    print("Iparams: ", so.iparam)
    print(so)
    print("name: ", so.name)
    so.print()
    assert so.iparam[pnames.SICONOS_IPARAM_MAX_ITER] == 1000
    assert so.solverId == solver_ids.SICONOS_FRICTION_3D_NSGS
    assert so.dparam[pnames.SICONOS_DPARAM_TOL] == 1.0e-4
    assert so.name == "FC3D_NSGS"
    so2 = so.get_internal_solver(0)
    print(f"Internal solver name: {so2.name} - Internal solver id: {so2.solverId}")
    assert so2.name == "ONECONE_NSN_GP_HYBRID"
    so.update_internal(0, 553)
    print("update ...")
    so2 = so.get_internal_solver(0)
    print(f"Internal solver name: {so2.name} - Internal solver id: {so2.solverId}")
    assert so2.name == "OC_PROJ_LI"


if __name__ == "__main__":
    test_solver_options_create()
    print("All tests passed!")
