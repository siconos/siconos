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

import siconos.mechanics.joints as mj
import siconos.modeling as sm
import numpy as np


def test_joints():
    j1 = mj.CylindricalJointR()
    j2 = mj.CylindricalJointR()
    j3 = mj.PrismaticJointR(np.ones(3), True)
    j4 = mj.PrismaticJointR(np.ones(3), True)
    assert isinstance(j1, mj.NewtonEulerJointR)
    cj = mj.CouplerJointR(j1, 0, j2, 1, 0.5, np.ones(7), 1, np.ones(7), 1)
    ck = mj.CouplerJointR(j3, 0, j4, 1, 0.5)
    q0 = np.ones(7, dtype=np.float64)
    twist0 = np.ones(6, dtype=np.float64)
    mass = 2.0
    inertia = np.zeros((3, 3), dtype=np.float64, order="F")
    inertia[0, 0] = inertia[1, 1] = inertia[2, 2] = 1.0
    neds = sm.NewtonEulerDS(q0, twist0, mass, inertia)

    q1 = neds.q()
    q2 = None
    ck.setBasePositions(q1, q2)
    neds2 = sm.NewtonEulerDS(q0, twist0, mass, inertia)
    cj.setBasePositions(q1, neds2.q())
