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
"""Python tests for Lagrangian dynamical systems classes and functions
   from kernel/modelingtools

"""

import cppimport.import_hook
import addons.computeLDS

import numpy as np
import siconos.modeling as sm
import scipy.sparse as sp
from scipy.sparse import csc_array

rng = np.random.default_rng(seed=42)


# -- Dense storage Lagrangian DS --
def test_LagrangianDS_alias():
    print("start test_LagrangianDS")
    q0 = np.array([1, 0, 0], dtype=np.float64)
    v0 = np.zeros_like(q0)
    ball = sm.LagrangianDS(q0, v0, sm.alias_t)

    assert np.allclose(q0, ball.q0)
    assert np.allclose(v0, ball.velocity0)

    q0[...] += 3
    v0[...] += 1.2

    assert np.allclose(q0, ball.q0)
    assert np.allclose(v0, ball.velocity0)

    fext = np.zeros_like(q0)
    fext[:] = 122
    ball.setConstantFext(fext, sm.alias_t)
    assert np.allclose(ball.fext[:], fext[:])
    fext[1] = 12
    assert np.allclose(ball.fext[:], fext[:])

    fext = 453
    assert ball.fext[0] == ball.fext[2] == 122
    assert ball.fext[1] == 12
    print("end test_LagrangianDS")


def test_LagrangianDS_copy():
    print("start test_LagrangianDS")
    q0 = np.array([1, 0, 0], dtype=np.float64)
    v0 = np.zeros_like(q0)
    ball = sm.LagrangianDS(q0, v0, sm.copy_t)

    assert np.allclose(q0, ball.q0)
    assert np.allclose(v0, ball.velocity0)

    q0[...] += 3
    v0[...] += 1.2

    assert not np.allclose(q0, ball.q0)
    assert not np.allclose(v0, ball.velocity0)

    fext = np.zeros_like(q0)
    fext[:] = 122
    ball.setConstantFext(fext, sm.alias_t)
    assert np.allclose(ball.fext[:], fext[:])
    fext[1] = 12
    assert np.allclose(ball.fext[:], fext[:])

    fext = 453
    assert ball.fext[0] == ball.fext[2] == 122
    assert ball.fext[1] == 12
    print("end test_LagrangianDS")


def fint_func(vel, pos, time, result):
    result[0] = vel[1] * np.cos(time)
    result[:] = vel * time + pos


def mass_func(pos, result):
    for i in range(pos.size):
        result[:, i] = i * pos


def test_LagrangianDS_compute_mass():
    print("start test_LagrangianDS_compute_mass")
    ndof = 3
    q0 = np.zeros(ndof, dtype=np.float64)
    v0 = np.zeros_like(q0)
    q0[:] = np.arange(1, ndof + 1)
    v0[:] = np.arange(4, ndof + 4)
    ball = sm.LagrangianDS(q0, v0, sm.alias_t)
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


def test_LagrangianDS_compute():
    print("start test_LagrangianDS_compute")

    ndof = 3
    q0 = np.zeros(ndof, dtype=np.float64)
    v0 = np.zeros_like(q0)
    q0[:] = np.arange(1, ndof + 1)
    v0[:] = np.arange(4, ndof + 4)
    ball = sm.LagrangianDS(q0, v0, sm.alias_t)
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


def test_lagrangianLinearTIDS_alias():
    print("start test_lagrangianLinearTIDS_alias")
    q0 = np.array([1, 0, 0], dtype=np.float64)
    v0 = np.zeros_like(q0)
    mass = np.empty((3, 3), dtype=np.float64, order="F")
    stiff = np.empty((3, 3), dtype=np.float64, order="F")
    damp = np.empty((3, 3), dtype=np.float64, order="F")
    rng.random(mass.shape, out=mass)
    rng.random(mass.shape, out=stiff)
    rng.random(mass.shape, out=damp)
    ball = sm.LagrangianLinearTIDS(q0, v0, mass, sm.alias_t)
    assert np.allclose(ball.q0, q0)
    assert np.allclose(ball.velocity0, v0)
    assert np.allclose(ball.mass, mass)
    q0 *= 3.0
    v0 += 1
    mass *= 2.1
    assert np.allclose(ball.q0, q0)
    assert np.allclose(ball.velocity0, v0)
    assert np.allclose(ball.mass, mass)

    ball.setStiffnessMatrix(stiff, sm.alias_t)
    ball.setDampingMatrix(damp, sm.alias_t)

    assert np.allclose(ball.stiffnessMatrix, stiff)
    assert np.allclose(ball.dampingMatrix, damp)
    damp *= 1.2
    stiff += 14.0
    assert np.allclose(ball.stiffnessMatrix, stiff)
    assert np.allclose(ball.dampingMatrix, damp)

    fext = np.zeros_like(q0)
    fext[:] = 122
    ball.setConstantFext(fext, sm.alias_t)
    assert np.allclose(ball.fext[:], fext[:])
    fext[1] = 12
    assert np.allclose(ball.fext[:], fext[:])

    fext = 453
    assert ball.fext[0] == ball.fext[2] == 122
    assert ball.fext[1] == 12
    print("end test_lagrangianLinearTIDS_alias")


def test_lagrangianLinearTIDS_copy():
    print("start test_lagrangianLinearTIDS_copy")
    q0 = np.array([1, 0, 0], dtype=np.float64)
    v0 = np.zeros_like(q0)
    mass = np.empty((3, 3), dtype=np.float64, order="F")
    stiff = np.empty((3, 3), dtype=np.float64, order="F")
    damp = np.empty((3, 3), dtype=np.float64, order="F")
    rng.random(mass.shape, out=mass)
    rng.random(mass.shape, out=stiff)
    rng.random(mass.shape, out=damp)
    ball = sm.LagrangianLinearTIDS(q0, v0, mass, sm.copy_t)
    assert np.allclose(ball.q0, q0)
    assert np.allclose(ball.velocity0, v0)
    assert np.allclose(ball.mass, mass)
    q0 *= 3.0
    v0 += 1
    mass *= 2.1
    assert not np.allclose(ball.q0, q0)
    assert not np.allclose(ball.velocity0, v0)
    assert not np.allclose(ball.mass, mass)

    ball.setStiffnessMatrix(stiff, sm.copy_t)
    ball.setDampingMatrix(damp, sm.copy_t)

    assert np.allclose(ball.stiffnessMatrix, stiff)
    assert np.allclose(ball.dampingMatrix, damp)
    damp *= 1.2
    stiff += 14.0
    assert not np.allclose(ball.stiffnessMatrix, stiff)
    assert not np.allclose(ball.dampingMatrix, damp)

    fext = np.zeros_like(q0)
    fext[:] = 122
    ball.setConstantFext(fext, sm.copy_t)
    assert np.allclose(ball.fext[:], fext[:])
    fext[1] = 12
    assert not np.allclose(ball.fext[:], fext[:])
    print("end test_lagrangianLinearTIDS_copy")


def test_NewtonEulerDS():
    q0 = np.ones(7, dtype=np.float64)
    twist0 = np.ones(6, dtype=np.float64)
    mass = 2.0
    inertia = np.zeros((3, 3), dtype=np.float64, order="F")
    inertia[0, 0] = inertia[1, 1] = inertia[2, 2] = 1.0
    neds = sm.NewtonEulerDS(q0, twist0, mass, inertia, sm.alias_t)
    assert np.allclose(neds.q(), q0)
    assert np.allclose(neds.twist(), twist0)


# -- Sparse storage Lagrangian DS --
def call_ds_compute(dstype, ndof):
    print(" ====== Compute test ====== ")
    initial_position = rng.random(ndof)
    initial_velocity = rng.random(ndof)

    ds = dstype(initial_position, initial_velocity, sm.alias_t)

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
    ds = dstype(initial_position, initial_velocity, sm.alias_t)

    ds.setConstantMass(mass, sm.alias_t)
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
    ds = dstype(initial_position, initial_velocity, sm.alias_t)

    ds.setConstantMass(mass, sm.copy_t)
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


def test_lagrangianSparseLinearTIDS_copy():
    print("start test_lagrangianSparseLinearTIDS_copy")
    q0 = np.array([1, 0, 0], dtype=np.float64)
    v0 = np.zeros_like(q0)
    data_m = rng.random(3)
    mass = sp.diags_array(data_m, format="csc")
    data_k = rng.random(3)
    stiff = sp.diags_array(data_k, format="csc")
    data_c = rng.random(3)
    damp = sp.diags_array(data_c, format="csc")

    ball = sm.LagrangianSparseLinearTIDS(q0, v0, mass, sm.copy_t)

    ball.setStiffnessMatrix(stiff, sm.copy_t)
    ball.setDampingMatrix(damp, sm.copy_t)
    fext = np.zeros_like(q0)
    fext[:] = 122
    ball.setConstantFext(fext, sm.alias_t)
    assert np.allclose(ball.fext[:], fext[:])
    fext[1] = 12
    assert np.allclose(ball.fext[:], fext[:])

    fext = 453
    assert ball.fext[0] == ball.fext[2] == 122
    assert ball.fext[1] == 12

    assert np.allclose(mass.toarray(), ball.mass_view.toarray())
    assert np.allclose(stiff.toarray(), ball.stiffness_view.toarray())
    assert np.allclose(damp.toarray(), ball.damping_view.toarray())
    stiff[1, 1] = 128
    assert not np.allclose(stiff.toarray(), ball.stiffness_view.toarray())
    print("end test_lagrangianSparseLinearTIDS_copy")


def test_lagrangianSparseLinearTIDS_alias():
    print("start test_lagrangianSparseLinearTIDS_alias")
    q0 = np.array([1, 0, 0], dtype=np.float64)
    v0 = np.zeros_like(q0)
    data_m = rng.random(3)
    mass = sp.diags_array(data_m, format="csc")
    data_k = rng.random(3)
    stiff = sp.diags_array(data_k, format="csc")
    data_c = rng.random(3)
    damp = sp.diags_array(data_c, format="csc")
    ball = sm.LagrangianSparseLinearTIDS(q0, v0, sm.alias_t)
    ball.setConstantMass(mass, sm.alias_t)
    ball.setStiffnessMatrix(stiff, sm.alias_t)
    ball.setDampingMatrix(damp, sm.alias_t)
    fext = np.zeros_like(q0)
    fext[:] = 122
    ball.setConstantFext(fext, sm.alias_t)
    assert np.allclose(ball.fext[:], fext[:])
    fext[1] = 12
    assert np.allclose(ball.fext[:], fext[:])

    fext = 453
    assert ball.fext[0] == ball.fext[2] == 122
    assert ball.fext[1] == 12

    assert np.allclose(mass.toarray(), ball.mass_alias.toarray())
    assert np.allclose(stiff.toarray(), ball.stiffness_alias.toarray())
    assert np.allclose(damp.toarray(), ball.damping_alias.toarray())
    stiff[1, 1] = 128
    assert np.allclose(stiff.toarray(), ball.stiffness_alias.toarray())
    print("end test_lagrangianSparseLinearTIDS_alias")


if __name__ == "__main__":
    test_LagrangianDS_alias()
    test_LagrangianDS_copy()
    test_LagrangianDS_compute()
    test_LagrangianDS_compute_mass()
    test_lagrangianLinearTIDS_alias()
    test_lagrangianLinearTIDS_copy()
    test_LagrangianSparseDS_copy()
    test_LagrangianSparseDS_alias()
    test_LagrangianSparseDS_compute()
    test_lagrangianSparseLinearTIDS_copy()
    test_lagrangianSparseLinearTIDS_alias()()
    test_NewtonEulerDS()
