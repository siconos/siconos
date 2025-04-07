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
"""A few tests for classes and functions from kernel/modelingtools

"""

import numpy as np
import siconos.modeling as sm


def test_lagrangianDS():
    print("start test_lagrangianDS")
    q0 = np.array([1, 0, 0], dtype=np.float64)
    v0 = np.zeros_like(q0)
    ball = sm.LagrangianDS(q0, v0)

    fext = np.zeros_like(q0)
    fext[:] = 122
    ball.setConstantFext(fext)
    assert np.allclose(ball.fext()[:], fext[:])
    fext[1] = 12
    assert np.allclose(ball.fext()[:], fext[:])

    fext = 453
    assert ball.fext()[0] == ball.fext()[2] == 122
    assert ball.fext()[1] == 12
    print("end test_lagrangianDS")


def fint_func(vel, pos, time, result):
    result[0] = vel[1] * np.cos(time)
    result[:] = vel * time + pos


def mass_func(pos, result):
    for i in range(pos.size):
        result[:, i] = i * pos


def test_lagrangianDS_compute_mass():
    print("start test_lagrangianDS_compute_mass")
    ndof = 3
    q0 = np.zeros(ndof, dtype=np.float64)
    v0 = np.zeros_like(q0)
    q0[:] = np.arange(1, ndof + 1)
    v0[:] = np.arange(4, ndof + 4)
    ball = sm.LagrangianDS(q0, v0)
    ball.setComputeMassFunction(mass_func)

    q = ball.q()
    v = ball.velocity()
    p = ball.p(1)

    q[:] = 1.0
    v[:] = 2.0
    p[:] = 3.0

    ref = np.zeros((ndof, ndof), dtype=np.float64)
    ref[:] = 1.0
    mass_func(q, ref)

    ball.computeMass(q)
    assert np.allclose(ball.mass, ref)


def test_lagrangianDS_compute():
    print("start test_lagrangianDS_compute")

    ndof = 3
    q0 = np.zeros(ndof, dtype=np.float64)
    v0 = np.zeros_like(q0)
    q0[:] = np.arange(1, ndof + 1)
    v0[:] = np.arange(4, ndof + 4)
    ball = sm.LagrangianDS(q0, v0)
    ball.setComputeFintFunction(fint_func)

    q = ball.q()
    v = ball.velocity()
    p = ball.p(1)
    time = 1.0

    q[:] = 1.0
    v[:] = 2.0
    p[:] = 3.0
    ref = np.zeros(q.size)
    fint_func(v, q, time, ref)

    ball.computeFint(v, q, 1.0)
    assert np.allclose(ball.fint, ref)


def test_NewtonEulerDS():
    q0 = np.ones(7, dtype=np.float64)
    twist0 = np.ones(6, dtype=np.float64)
    mass = 2.0
    inertia = np.zeros((3, 3), dtype=np.float64, order="F")
    inertia[0, 0] = inertia[1, 1] = inertia[2, 2] = 1.0
    neds = sm.NewtonEulerDS(q0, twist0, mass, inertia)
    assert np.allclose(neds.q(), q0)
    assert np.allclose(neds.twist(), twist0)


if __name__ == "__main__":
    test_lagrangianDS()
    test_lagrangianDS_compute()
    test_lagrangianDS_compute_mass()
