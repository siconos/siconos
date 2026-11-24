# Siconos is a program dedicated to modeling, simulation and control
# of non smooth dynamical systems.
#
# Copyright INRIA.
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

"""energy_work (output_energy_work), fem_velocities_, fem_mass_ and
fem_dof_indices_ with the vnative backend.

The energy_work columns of a row are, for the step that ends at t: the
kinetic energy at t, the work of the forces h f_theta . v_theta (f the total
forces: fext, and - K q for a FEM body; theta of the integrator), and the
works of the contact impulses (1/2 and theta weighted). With M dv = h
f_theta + H^T lambda, the step closes

    K(k + 1) - K(k) = force_work - (theta - 1/2) h f_theta . dv + w_n + w_t

Run with pytest, nonos and siconos.io in the PYTHONPATH; the FEM tests take
the mesh mesh_N5.msh of this directory.
"""

import os

import h5py
import numpy as np
import pytest
import scipy.constants
import scipy.sparse as sp

import siconos.numerics as sn
from siconos.mechanics.collision.tools import Contactor
from siconos.io.mechanics_run import (
    MechanicsHdf5Runner,
    MechanicsHdf5Runner_run_options,
    RunnerConfig,
)
import nonos  # noqa: F401 registers the vnative backend

G = scipy.constants.g
H = 0.005
MESH = os.path.join(os.path.dirname(os.path.abspath(__file__)), "mesh_N5.msh")


def no_forces(body):
    pass


def run(filename, theta, T, **options):
    run_options = MechanicsHdf5Runner_run_options()
    run_options["t0"] = 0
    run_options["T"] = T
    run_options["h"] = H
    run_options["theta"] = theta
    so = sn.SolverOptions(sn.solver_ids.SICONOS_FRICTION_2D_NSGS)
    so.iparam[sn.params.SICONOS_IPARAM_MAX_ITER] = 1000
    so.dparam[sn.params.SICONOS_DPARAM_TOL] = 1e-12
    run_options["solver_options"] = so
    run_options["output_energy_work"] = True
    run_options["output_frequency"] = 1
    run_options["verbose"] = False
    run_options.update(options)
    with MechanicsHdf5Runner(
        mode="r+", config=RunnerConfig("vnative"), io_filename=filename
    ) as io:
        io.run(run_options)
    with h5py.File(filename, "r") as f:
        return {
            name: f["data"][name][()]
            for name in f["data"]
            if name != "log" and isinstance(f["data"][name], h5py.Dataset)
        }


def disk_scene(filename, velocity, ground=False, e=0.0):
    with MechanicsHdf5Runner(
        mode="w", config=RunnerConfig("vnative"), io_filename=filename
    ) as io:
        io.add_primitive_shape("Disk", "Disk", [0.1])
        if ground:
            io.add_primitive_shape("Ground", "Segment", (-10, 0, 10, 0))
            io.add_object("ground", [Contactor("Ground")], translation=[0, 0])
        io.add_Newton_impact_friction_nsl("contact", mu=0.3, e=e)
        io.add_object(
            "disk",
            [Contactor("Disk")],
            translation=[0, 0.2],
            orientation=[0],
            velocity=velocity,
            mass=2.0,
            inertia=0.01,
        )


@pytest.mark.parametrize("theta", [1.0, 0.50001])
def test_free_fall(tmp_path, theta):
    """K = 1/2 m (g t)^2; force_work - dK = (theta - 1/2) m g^2 h^2, the
    numerical dissipation of the theta-method, not a rounding; no contact."""
    filename = str(tmp_path / "free_fall.hdf5")
    disk_scene(filename, [0, 0, 0])
    data = run(filename, theta, T=0.1)
    ew = data["energy_work"]
    t, kinetic, force_work = ew[:, 0], ew[:, 1], ew[:, 2]
    m = 2.0
    assert len(t) == 20
    assert kinetic == pytest.approx(0.5 * m * (G * t) ** 2, rel=1e-6)
    dk = np.diff(np.concatenate([[0.0], kinetic]))
    assert force_work - dk == pytest.approx(
        np.full(len(t), (theta - 0.5) * m * G**2 * H**2), abs=2e-6
    )
    assert np.all(ew[:, 3:] == 0)


def test_rotating_disk(tmp_path):
    """No force: K = 1/2 m v^2 + 1/2 I omega^2, constant, no work."""
    filename = str(tmp_path / "rotating.hdf5")
    disk_scene(filename, [0.3, 0, 5.0])
    data = run(filename, 0.50001, T=0.05, set_external_forces=no_forces)
    ew = data["energy_work"]
    assert ew[:, 1] == pytest.approx(0.5 * 2.0 * 0.3**2 + 0.5 * 0.01 * 5.0**2, rel=1e-6)
    assert np.all(ew[:, 2] == 0)


@pytest.mark.parametrize("e", [0.0, 0.5])
@pytest.mark.parametrize("theta", [0.5, 0.50001])
def test_bouncing_disk_balance(tmp_path, e, theta):
    """A disk dropped on a segment: each step closes the energy balance with
    the contact works (the normal of a segment does not turn)."""
    filename = str(tmp_path / "bouncing.hdf5")
    disk_scene(filename, [0, 0, 0], ground=True, e=e)
    data = run(filename, theta, T=0.6)
    ew = data["energy_work"]
    kinetic, force_work, w_n, w_t = ew[:, 1], ew[:, 2], ew[:, 3], ew[:, 4]
    # the velocity of the disk at each output step (one body)
    vy = data["velocities"][:, 3]
    v = np.concatenate([[0.0], vy[1 : len(kinetic) + 1]])
    dv = np.diff(v)
    dissipation = (theta - 0.5) * H * (-2.0 * G) * dv
    dk = np.diff(np.concatenate([[0.0], kinetic]))
    assert np.abs(w_n).max() > 0.1  # an impact happened
    residual = dk - force_work + dissipation - w_n - w_t
    assert np.abs(residual).max() < 1e-5 * np.abs(kinetic).max()


def applied_force_node():
    """the coordinates of the node of the "Applied Force" tag (3) of the mesh"""
    lines = open(MESH).read().split("\n")
    nodes = lines.index("$Nodes")
    elements = lines.index("$Elements")
    for line in lines[elements + 2 :]:
        fields = line.split()
        if len(fields) > 4 and fields[1] == "15" and fields[3] == "3":
            node = int(fields[-1])
            x, y = map(float, lines[nodes + 1 + node].split()[1:3])
            return x, y
    raise RuntimeError("no Applied Force node")


def fem_scene(filename, rho, force, boundary_conditions):
    with MechanicsHdf5Runner(
        mode="w", config=RunnerConfig("vnative"), io_filename=filename
    ) as io:
        io.add_shape_data_from_file("Square", MESH)
        io.add_object(
            "square",
            [Contactor("Square")],
            material=(rho, 1e7, 0.3),
            nodal_forces=force,
            boundary_conditions=boundary_conditions,
        )


def read_csc(group):
    shape = (int(group.attrs["rows"]), int(group.attrs["cols"]))
    return sp.csc_array(
        (group["values"][:], group["inner"][:], group["outer"][:]), shape
    )


def fem_energies(tmp_path, theta, boundary_conditions):
    """A square, a constant force on a top corner from rest, no contact:
    the file read back (M, K, the dof map, u and v in spatial order), and
    the energy E = K + 1/2 u^T K u - f^T u at each step."""
    filename = str(tmp_path / "fem.hdf5")
    rho, force = 2500.0, np.array([1e4, -2e4])
    fem_scene(filename, rho, force, boundary_conditions)
    data = run(filename, theta, T=0.1, set_external_forces=no_forces)
    with h5py.File(filename, "r") as f:
        d = f["data"]
        (ds_id,) = [
            int(k.removeprefix("fem_mass_")) for k in d if k.startswith("fem_mass_")
        ]
        M = read_csc(d[f"fem_mass_{ds_id}"])
        K = read_csc(d[f"fem_stiffness_{ds_id}"])
        assert d[f"fem_mass_{ds_id}"].attrs["mass"] == "consistent"
        dofs = d[f"fem_dof_indices_{ds_id}"][:]
        coords = d[f"fem_coords_{ds_id}"][:]
        u_rows = d[f"fem_displacements_{ds_id}"][:]
        v_rows = d[f"fem_velocities_{ds_id}"][:]
        assert d[f"fem_velocities_{ds_id}"].dtype == np.float64
        assert dict(d[f"fem_velocities_{ds_id}"].attrs) == dict(
            d[f"fem_displacements_{ds_id}"].attrs
        )

    n = M.shape[0]
    # M: symmetric, the mass of the unit square (thickness 1) in each
    # direction, consistent (not diagonal)
    assert abs(M - M.T).max() < 1e-12 * abs(M).max()
    assert M.sum() == pytest.approx(2 * rho, rel=1e-12)
    assert M.nnz > n
    assert sorted(dofs) == list(range(n))

    # the FEM tables start at the first step: the initial state, at rest and
    # undeformed, is row 0 here
    def full(rows):
        x = np.zeros((rows.shape[0] + 1, n))
        x[1:, dofs] = rows[:, 2:]
        return x

    u, v = full(u_rows), full(v_rows)
    assert u_rows[:, 0] == pytest.approx(v_rows[:, 0], rel=1e-6)
    # the positions follow the theta-method: du / h = theta v(k+1) + (1-theta) v(k)
    # (fem_displacements_ is float32 in the file, fem_velocities_ float64)
    eps32 = np.finfo(np.float32).eps
    assert np.diff(u, axis=0) / H == pytest.approx(
        theta * v[1:] + (1 - theta) * v[:-1], rel=0, abs=2 * eps32 * np.abs(u).max() / H
    )

    ew = data["energy_work"]
    kinetic = np.array([0.5 * x @ (M @ x) for x in v])
    # 1/2 v^T M v from the file is the kinetic column
    assert ew[:, 1] == pytest.approx(kinetic[1:], rel=1e-6, abs=1e-6 * kinetic.max())

    # f: the force on the applied force node
    (i,) = np.flatnonzero(
        np.all(np.isclose(coords[:, :2], applied_force_node()), axis=1)
    )
    f = np.zeros(n)
    f[dofs[2 * i : 2 * i + 2]] = force
    energy = kinetic + 0.5 * np.einsum("ij,ij->i", u, (K @ u.T).T) - u @ f
    # the dissipation of the theta-method, (theta - 1/2) h f_theta . dv
    u_theta = theta * u[1:] + (1 - theta) * u[:-1]
    f_theta = f - (K @ u_theta.T).T
    dissipation = (theta - 0.5) * H * np.einsum("ij,ij->i", f_theta, np.diff(v, axis=0))
    return energy, ew[:, 2], np.diff(kinetic), dissipation, np.abs(u @ f).max()


@pytest.mark.parametrize("theta", [0.5, 0.6])
def test_fem_free_oscillation(tmp_path, theta):
    """A free square (no Dirichlet condition) under a constant nodal force:
    K + 1/2 u^T K u - f^T u constant at theta = 1/2, decaying at theta >
    1/2; force_work (the elastic work included) closes dK per step."""
    energy, force_work, dk, dissipation, scale = fem_energies(tmp_path, theta, [])
    assert scale > 0
    assert force_work - dk == pytest.approx(dissipation, abs=1e-6 * scale)
    if theta == 0.5:
        assert np.abs(energy - energy[0]).max() < 1e-6 * scale
    else:
        assert np.all(np.diff(energy) < 1e-7 * scale)
        assert energy[-1] < energy[0] - 1e-3 * scale


@pytest.mark.parametrize("theta", [0.5, 0.6])
def test_fem_clamped_energy(tmp_path, theta):
    """The same square clamped at its base (the Dirichlet nodes of the mesh):
    the free dofs are solved with the operator of the free dofs, and the
    balance closes as for the free square (the blocked dofs do no work)."""
    energy, force_work, dk, dissipation, scale = fem_energies(tmp_path, theta, None)
    assert scale > 0
    assert force_work - dk == pytest.approx(dissipation, abs=1e-6 * scale)
    if theta == 0.5:
        assert np.abs(energy - energy[0]).max() < 1e-6 * scale
    else:
        assert np.all(np.diff(energy) < 1e-7 * scale)
        assert energy[-1] < energy[0] - 1e-3 * scale


def test_fem_restart(tmp_path):
    """A cold restart (run again to a later T) appends to energy_work,
    fem_velocities_ and fem_displacements_, keeps one fem_mass_,
    fem_stiffness_, fem_dof_indices_ (their metadata checked), and the FEM
    body goes on from its last output state: the same trajectory as one run,
    up to the float32 of fem_displacements_ the restart reads."""
    force = np.array([1e4, -2e4])
    tables = {}
    for name, ends in (("one_run", [0.1]), ("restarted", [0.05, 0.1])):
        filename = str(tmp_path / f"{name}.hdf5")
        fem_scene(filename, 2500.0, force, None)
        for T in ends:
            run(filename, 0.5, T=T, set_external_forces=no_forces)
        with h5py.File(filename, "r") as f:
            d = f["data"]
            groups = sorted(
                k for k in d if k.startswith(("fem_mass_", "fem_stiffness_", "fem_dof"))
            )
            assert [g.rsplit("_", 1)[0] for g in groups] == [
                "fem_dof_indices",
                "fem_mass",
                "fem_stiffness",
            ]
            (ds_id,) = {g.rsplit("_", 1)[1] for g in groups}
            tables[name] = {
                what: d[f"{what}_{ds_id}" if what != "energy_work" else what][:]
                for what in ("fem_velocities", "fem_displacements", "energy_work")
            }
    for what, one_run in tables["one_run"].items():
        restarted = tables["restarted"][what]
        assert restarted.shape == one_run.shape == (20, one_run.shape[1])
        assert restarted[:, 0] == pytest.approx(H * np.arange(1, 21), rel=1e-6)
        scale = np.abs(one_run[:, 1:]).max()
        assert np.abs(restarted[:, 1:] - one_run[:, 1:]).max() < 1e-5 * scale
