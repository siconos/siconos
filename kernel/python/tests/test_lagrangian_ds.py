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
"""Python tests for Lagrangian dynamical systems classes and functions from kernel/modelingtools

"""

import cppimport.import_hook
import addons.computeLDS

import numpy as np
import siconos.modeling as sm
import scipy.sparse as sp
from scipy.sparse import csc_array

rng = np.random.default_rng(seed=42)


# -- Dense storage Lagrangian DS --
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


# -- Sparse storage Lagrangian DS --
def call_ds_compute(dstype, ndof):
    print(" ====== Compute test ====== ")
    initial_position = rng.random(ndof)
    initial_velocity = rng.random(ndof)

    ds = dstype(initial_position, initial_velocity)

    def ref_func(q):
        diag_coeffs = np.array(1.0 + 0.1 * q)
        M = csc_array(sp.diags_array(diag_coeffs))
        return M

    # Compute a reference mass
    vec = np.array([1.0, 22.0, 103.0], dtype=np.float64)
    Mref = ref_func(vec)
    # Initialize DS mass
    ds.setComputeMassFunction(addons.computeLDS.initializeSparse)
    ds.computeMass(vec)
    print("Mref:\n", Mref.toarray())
    print("Initial mass value\n", ds.mass_view.toarray())
    # print(ds)

    assert np.allclose(Mref.toarray(), ds.mass_view.toarray())

    vec = rng.random(ds.dimension)
    ref = Mref @ vec
    result = ds.mass_view @ vec
    assert np.allclose(ref, result)

    M = ds.mass_view
    print("M\n", M.toarray())
    print("mass\n", ds.mass_view.toarray())

    assert not ds.mass_view.data.flags.writeable

    vec[0] *= 12
    Mref2 = Mref
    # Compute Mass
    ds.setComputeMassFunction(addons.computeLDS.computeMassSparse)
    ds.computeMass(vec)
    Mref *= 1000
    assert np.allclose(Mref.toarray(), ds.mass_view.toarray())
    assert np.allclose(Mref.toarray(), Mref2.toarray())

    print("mview end", ds.mass_view)


def call_ds_alias(dstype, ndof):
    print(" ====== Alias test ====== ")
    initial_position = rng.random(ndof)
    initial_velocity = rng.random(ndof)

    data = rng.random(ndof)
    mass = sp.diags_array(data, format="csc")
    ds = dstype(initial_position, initial_velocity)

    ds.setConstantMassAlias(mass)
    try:
        M = ds.mass_view
        assert False
    except:
        print("ok")
        pass

    M = ds.mass_alias
    mass[0, 0] = -99.0

    assert np.allclose(mass.toarray(), ds.mass_alias.toarray())
    assert np.allclose(M.toarray(), ds.mass_alias.toarray())
    vec = rng.random(ds.dimension)
    ref = mass @ vec
    result = ds.mass_alias @ vec
    assert np.allclose(ref, result)


def call_ds_copy(dstype, ndof):
    print(" ====== Copy test ====== ")
    initial_position = rng.random(ndof)
    initial_velocity = rng.random(ndof)

    data = rng.random(ndof)
    mass = sp.diags_array(data, format="csc")
    ds = dstype(initial_position, initial_velocity)

    ds.setConstantMassCopy(mass)
    try:
        M = ds.mass_alias
        assert False

    except:
        print("ok")
        pass

    M = ds.mass_view
    vec = rng.random(ds.dimension)
    res = M @ vec
    ref = mass @ vec
    assert np.allclose(ref, res)
    mass[0, 0] = -99.0
    assert not np.allclose(mass.toarray(), ds.mass_view.toarray())


def test_LagrangianSparseDS_copy():
    ndof = 3
    call_ds_copy(sm.LagrangianSparseDS, ndof)


def test_LagrangianSparseDS_alias():
    ndof = 3
    call_ds_alias(sm.LagrangianSparseDS, ndof)


def test_LagrangianSparseDS_compute():
    ndof = 3
    call_ds_compute(sm.LagrangianSparseDS, ndof)


if __name__ == "__main__":
    test_lagrangianDS()
    test_lagrangianDS_compute()
    test_lagrangianDS_compute_mass()
    test_LagrangianSparseDS_copy()
    test_LagrangianSparseDS_alias()
    test_LagrangianSparseDS_compute()
