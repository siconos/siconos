# Siconos is a program dedicated to modeling, simulation and control
# of non smooth dynamical systems.
#
# SPDX-FileCopyrightText: INRIA
# SPDX-License-Identifier: Apache-2.0
#
# Check: a small pile of disks under the quasi inelastic shock law, against
# the Newton law (e = 0).  A 3 x 3 grid of disks falls under gravity on a
# floor and settles.  The law settles the pile without the ejection chattering
# of the velocity Signorini condition, and the contacts that appear because a
# gap closes are kept (the active contact count is stable once settled).
#
# usage: python3 quasi_inelastic_pile.py
# exit status 0 if the quasi inelastic pile settles and does not eject

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

R, m = 0.1, 1.0
inertia = 0.5 * m * R * R
N = 3  # 3 x 3
T = 0.6


def run(osi_name, filename):
    config = RunnerConfig("vnative")
    spacing = 2.05 * R
    pos = []
    for j in range(N):
        for i in range(N):
            pos.append([(i - (N - 1) / 2) * spacing, R + (j + 1) * spacing * 1.0])
    pos = np.array(pos)
    with MechanicsHdf5Runner(config=config, io_filename=filename) as io:
        io.add_primitive_shape("Disk", "Disk", [R])
        io.add_primitive_shape("Ground", "Segment", (-10, 0, 10, 0))
        io.add_objects(
            "disks",
            [Contactor("Disk")],
            translations=pos,
            orientations=np.zeros(len(pos)),
            velocities=np.zeros((len(pos), 3)),
            mass=m,
            inertia=inertia,
        )
        io.add_object("ground", [Contactor("Ground")], translation=[0, 0])
        io.add_Newton_impact_friction_nsl("contact", mu=0.5, e=0.0)

    options = sn.SolverOptions(sn.solver_ids.SICONOS_FRICTION_2D_NSGS)
    options.iparam[sn.params.SICONOS_IPARAM_MAX_ITER] = 1000
    options.dparam[sn.params.SICONOS_DPARAM_TOL] = 1e-10

    ro = MechanicsHdf5Runner_run_options()
    ro["t0"] = 0.0
    ro["T"] = T
    ro["h"] = 2e-3
    ro["theta"] = 0.5
    ro["solver_options"] = options
    if osi_name == "qi":
        ro["osi"] = nonos.bridge.MoreauJeanQuasiInelasticOSI
    with MechanicsHdf5Runner(mode="r+", config=config, io_filename=filename) as io:
        io.run(ro)

    with h5py.File(filename, "r") as f:
        vel = f["data/velocities"][:]
        cf = f["data/cf"][:]
    return vel, cf


def active_per_step(cf):
    if cf.shape[0] == 0:
        return np.zeros(0, dtype=int)
    times, counts = np.unique(cf[:, 0], return_counts=True)
    return counts


if len(sys.argv) == 3:
    vel, cf = run(sys.argv[1], sys.argv[2])
    np.savez(sys.argv[2] + ".npz", vel=vel, cf=cf)
    sys.exit(0)

tmpdir = tempfile.mkdtemp()
data = {}
for name in ("qi", "newton"):
    out = os.path.join(tmpdir, name)
    subprocess.run([sys.executable, __file__, name, out], check=True)
    data[name] = np.load(out + ".npz")

failed = False


def settled_stats(name):
    d = data[name]
    vel, cf = d["vel"], d["cf"]
    t = vel[:, 0]
    late = t > (T - 0.1)
    speed = np.hypot(vel[late, 2], vel[late, 3])
    up = vel[late, 3].max()
    counts = active_per_step(cf)
    chatter = int(np.abs(np.diff(counts)).sum()) if counts.size > 1 else 0
    print(
        f"{name}: max |v| late={speed.max():.3e} max upward vy={up:.3e} "
        f"active contacts {counts.min() if counts.size else 0}.."
        f"{counts.max() if counts.size else 0} chatter={chatter}"
    )
    return speed.max(), up, counts, chatter


sq, uq, cq, chq = settled_stats("qi")
sn_, un, cn, chn = settled_stats("newton")

# the quasi inelastic pile settles and does not eject
if sq > 1e-2:
    print(f"FAILED: the quasi inelastic pile does not settle (max |v| {sq:.3e})")
    failed = True
if uq > 1e-2:
    print(f"FAILED: the quasi inelastic pile ejects (upward vy {uq:.3e})")
    failed = True
# once settled the active contact count is stable (the contacts are kept, not
# recreated step after step)
if cq.size > 5:
    late_c = cq[-5:]
    if late_c.max() - late_c.min() > 1:
        print(f"FAILED: the active contact count chatters once settled ({late_c})")
        failed = True

sys.exit(1 if failed else 0)
