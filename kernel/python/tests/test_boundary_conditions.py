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
"""Python tests for BoundaryConditions from kernel/modeling

"""

import numpy as np
import siconos.modeling as sm

rng = np.random.default_rng(seed=42)


def test_bc_fixed():
    bc = sm.BoundaryCondition([1, 2, 4])

    assert bc.size == 3
    assert np.allclose(bc.prescribedVelocity, np.zeros(3))
    assert np.allclose(bc.velocityIndices, [1, 2, 4])


def test_bc_with_values():
    velo = np.asarray([1.0, 2.0, 3.0], dtype=np.float64)
    bc = sm.BoundaryCondition([1, 2, 4], velo)
    assert bc.size == 3
    assert np.allclose(bc.prescribedVelocity, velo)
    assert np.allclose(bc.velocityIndices, [1, 2, 4])


def test_harmonic_bc_scalar():
    a = 1.0
    b = 2.0
    omega = 3.0
    phi = 4.0
    bc = sm.HarmonicBC([1, 2, 4], a, b, omega, phi)
    time = 2.0
    bc.computePrescribedVelocity(time)

    velo = np.zeros(3)
    velo[...] = a + b * np.cos(omega * time + phi)

    assert bc.size == 3
    assert np.allclose(bc.prescribedVelocity, velo)
    assert np.allclose(bc.velocityIndices, [1, 2, 4])


def test_harmonic_bc_vector():
    size = 4
    a = np.empty(size, dtype=np.float64)
    rng.random((size,), out=a)
    b = np.empty(size, dtype=np.float64)
    rng.random((size,), out=b)
    omega = np.empty(size, dtype=np.float64)
    rng.random((size,), out=omega)
    phi = np.empty(size, dtype=np.float64)
    rng.random((size,), out=phi)
    bc = sm.HarmonicBC([1, 2, 4, 8], a, b, omega, phi)
    time = 2.0
    bc.computePrescribedVelocity(time)

    velo = a + b * np.cos(omega * time + phi)

    assert bc.size == 4
    assert np.allclose(bc.prescribedVelocity, velo)
    assert np.allclose(bc.velocityIndices, [1, 2, 4, 8])


if __name__ == "__main__":
    test_bc_fixed()
    test_bc_with_values()
    test_harmonic_bc_scalar()
    test_harmonic_bc_vector()
