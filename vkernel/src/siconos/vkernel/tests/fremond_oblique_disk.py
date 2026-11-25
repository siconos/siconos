# SPDX-FileCopyrightText: INRIA
#
# SPDX-License-Identifier: Apache-2.0
#
# Siconos is a program dedicated to modeling, simulation and control
# of non smooth dynamical systems.

# Check: the Frémond impact law and the inertia of the disks, on the vnative
# backend, through the runner.
#
# A disk (R = 0.1, m = 1, I = 2/5 m R^2) hits a ground segment obliquely,
# v = (1, -1), no spin, e = 1, mu = 1, theta = 1/2: the first impact sticks,
# whatever the normal velocity (gravity is on). With W_TT = 7 / (2 m):
#
#   Newton:  u_T(k + 1) = 0                         vx+ = 5/7
#   Frémond: u_T(k + 1/2) = 0, the slip reverses    vx+ = 3/7, KE conserved
#
# (Acary & Collins-Craft, JTCAM 2025). 5/7 needs the given inertia: the one of
# a uniform disk, m R^2 / 2, gives 2/3. The cases:
#
#   newton         add_object, the Newton law, Moreau-Jean           5/7
#   scheme         the same, --scheme MoreauJeanFremond (simoptions)  3/7
#   fremond        add_object, the Frémond law                        3/7
#   fremond-at     the Frémond law on Anitescu-Tasora: an error
#   objects        add_objects, a shared inertia                      5/7
#   objects-each   add_objects, one inertia per grain                 5/7
#
# Each case runs in its own process: a second simulation in the same
# process crashes the vnative backend.
#
# usage: python3 fremond_oblique_disk.py
# exit status 0 if every case gives its value

import os
import subprocess
import sys
import tempfile

import h5py
import numpy as np

import siconos.numerics as sn
from siconos.io.mechanics_run import (
    MechanicsHdf5Runner,
    MechanicsHdf5Runner_run_options,
    RunnerConfig,
)
from siconos.mechanics.collision.tools import Contactor

import nonos  # noqa: F401 registers the vnative backend
import nonos.bridge

radius, m = 0.1, 1.0
inertia = 2 / 5 * m * radius**2
v0 = [1.0, -1.0, 0.0]
y0 = radius + 0.0105


def run(case, filename):
    config = RunnerConfig("vnative")
    with MechanicsHdf5Runner(config=config, io_filename=filename) as io:
        io.add_primitive_shape("Disk", "Disk", [radius])
        io.add_primitive_shape("Ground", "Segment", (-10, 0, 10, 0))
        if case.startswith("objects"):
            # two disks 1 apart, the same impact
            io.add_objects(
                "disks",
                Contactor("Disk"),
                [[0, y0], [1, y0]],
                velocities=[v0, v0],
                mass=m,
                inertia=inertia if case == "objects" else [inertia, inertia],
            )
        else:
            io.add_object(
                "disk",
                [Contactor("Disk")],
                translation=[0, y0],
                orientation=[0],
                velocity=v0,
                mass=m,
                inertia=inertia,
            )
        io.add_object("ground", [Contactor("Ground")], translation=[0, 0])
        if case.startswith("fremond"):
            io.add_Fremond_impact_friction_nsl("contact", mu=1.0, e=1.0)
        else:
            io.add_Newton_impact_friction_nsl("contact", mu=1.0, e=1.0)

    if case == "scheme":
        # what the siconos command does with --scheme MoreauJeanFremond
        from nonos import simoptions

        simoptions.apply({"scheme": "MoreauJeanFremond"})

    options = sn.SolverOptions(sn.solver_ids.SICONOS_FRICTION_2D_NSGS)
    options.iparam[sn.params.SICONOS_IPARAM_MAX_ITER] = 1000
    options.dparam[sn.params.SICONOS_DPARAM_TOL] = 1e-14

    run_options = MechanicsHdf5Runner_run_options()
    run_options["t0"] = 0
    run_options["T"] = 0.02
    run_options["h"] = 1e-3
    run_options["theta"] = 0.5
    run_options["solver_options"] = options
    if case == "fremond-at":
        run_options["osi"] = nonos.bridge.AnitescuTasoraOSI
    with MechanicsHdf5Runner(mode="r+", config=config, io_filename=filename) as io:
        io.run(run_options)

    # rows (t, id, vx, vy, vz, wx, wy, wz) of data/velocities
    with h5py.File(filename, "r") as f:
        return f["data/velocities"][:]


if len(sys.argv) == 3:
    case, out = sys.argv[1:]
    try:
        np.save(out, run(case, out + ".hdf5"))
    except ValueError as e:
        print(f"ValueError: {e}")
        sys.exit(3)
    sys.exit(0)

expected = {
    "newton": 5 / 7,
    "scheme": 3 / 7,
    "fremond": 3 / 7,
    "objects": 5 / 7,
    "objects-each": 5 / 7,
}
# the velocities are float32 in the hdf5 file
tol = 1e-6
failed = False
tmpdir = tempfile.mkdtemp()
for case, vx_ref in expected.items():
    out = os.path.join(tmpdir, case + ".npy")
    subprocess.run(
        [sys.executable, __file__, case, out], check=True, stdout=subprocess.DEVNULL
    )
    v = np.load(out)
    for ds in np.unique(v[:, 1]):
        rows = v[v[:, 1] == ds]
        before = np.argmax(rows[:, 3] > 0) - 1  # the last step before the impact
        vx, vy, w = rows[before + 1, [2, 3, 7]]
        ke = 0.5 * m * (vx**2 + vy**2) + 0.5 * inertia * w**2
        ke_before = 0.5 * m * (rows[before, 2] ** 2 + rows[before, 3] ** 2)
        print(
            f"{case} ds {int(ds)}: vx+ = {vx:.9f} (expected {vx_ref:.9f}), "
            f"KE {ke_before:.9f} -> {ke:.9f}"
        )
        if abs(vx - vx_ref) > tol or abs(vy + rows[before, 3]) > tol:
            print(f"FAILED: {case}, vx+ or vy+")
            failed = True
        if vx_ref == 3 / 7 and abs(ke - ke_before) > 1e-5:
            print(f"FAILED: {case}, the kinetic energy is not conserved")
            failed = True

# the Frémond law on Anitescu-Tasora: the scheme does not implement it
out = os.path.join(tmpdir, "fremond-at.npy")
r = subprocess.run(
    [sys.executable, __file__, "fremond-at", out], capture_output=True, text=True
)
print(f"fremond-at: exit status {r.returncode}, {r.stdout.strip()}")
if r.returncode != 3 or "needs a MoreauJeanOSI" not in r.stdout:
    print("FAILED: fremond-at should be refused")
    failed = True

sys.exit(1 if failed else 0)
