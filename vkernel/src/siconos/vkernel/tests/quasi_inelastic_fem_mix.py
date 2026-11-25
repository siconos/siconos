# Siconos is a program dedicated to modeling, simulation and control
# of non smooth dynamical systems.
#
# SPDX-FileCopyrightText: INRIA
# SPDX-License-Identifier: Apache-2.0
#
# Check: a world that mixes a FEM mesh and a rigid disk, under the quasi
# inelastic shock law.  The law applies to the rigid-rigid contacts only (Jean
# sec 4.2); the FEM contact keeps the Newton impact law, and the setup warns.
# The run must complete and print the warning once.
#
# usage: python3 quasi_inelastic_fem_mix.py
# exit status 0 if the warning is printed and the run completes

import os
import subprocess
import sys
import tempfile

import siconos.numerics as sn
from siconos.io.mechanics_run import (
    MechanicsHdf5Runner,
    MechanicsHdf5Runner_run_options,
    RunnerConfig,
)
from siconos.mechanics.collision.tools import Contactor

import nonos  # noqa: F401 registers the vnative backend
import nonos.bridge

WARNING = "rigid-rigid contacts only"


def run(mesh, filename):
    config = RunnerConfig("vnative")
    with MechanicsHdf5Runner(config=config, io_filename=filename) as io:
        io.add_primitive_shape("DiskR", "Disk", [0.05])
        io.add_shape_data_from_file("Square", mesh)
        io.add_object("square", [Contactor("Square")], material=(2500, 1e11, 0.3))
        io.add_object(
            "disk",
            [Contactor("DiskR")],
            translation=[0.5, 1.0 + 0.05 + 0.01],
            orientation=[0],
            velocity=[0.0, -1.0, 0.0],
            mass=1,
            inertia=0.5 * 0.05 * 0.05,
        )
        io.add_Newton_impact_friction_nsl("contact", mu=0.5, e=0.0)

    options = sn.SolverOptions(sn.solver_ids.SICONOS_FRICTION_2D_NSGS)
    options.iparam[sn.params.SICONOS_IPARAM_MAX_ITER] = 100
    options.dparam[sn.params.SICONOS_DPARAM_TOL] = 1e-2
    options.iparam[sn.params.SICONOS_FRICTION_2D_SCALING_NORM] = (
        sn.params.SICONOS_FRICTION_2D_SCALING_NO
    )

    ro = MechanicsHdf5Runner_run_options()
    ro["t0"] = 0.0
    ro["T"] = 0.2
    ro["h"] = 5e-3
    ro["theta"] = 0.5
    ro["solver_options"] = options
    ro["osi"] = nonos.bridge.MoreauJeanQuasiInelasticOSI
    with MechanicsHdf5Runner(mode="r+", config=config, io_filename=filename) as io:
        io.run(ro)


if len(sys.argv) == 3:
    # the child: run the scene, the warning goes to stderr
    run(sys.argv[1], sys.argv[2])
    sys.exit(0)

tmpdir = tempfile.mkdtemp()
mesh = os.path.join(tmpdir, "square.msh")
subprocess.run(
    [
        sys.executable,
        os.path.join(os.path.dirname(__file__), "square_mesh.py"),
        "6",
        mesh,
    ],
    check=True,
)
r = subprocess.run(
    [sys.executable, __file__, mesh, os.path.join(tmpdir, "fem.hdf5")],
    capture_output=True,
    text=True,
)
print(f"exit {r.returncode}; warning present: {WARNING in r.stderr}")
if r.returncode != 0:
    print(r.stderr[-2000:])
    sys.exit(1)
if WARNING not in r.stderr:
    print("FAILED: the FEM-mix warning was not printed")
    sys.exit(1)
sys.exit(0)
