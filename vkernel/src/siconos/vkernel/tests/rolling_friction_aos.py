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

# Check: rolling resistance on the vnative backend against the AOS
# kernel. The reference is disk_line_rolling_aos (Lagrangian2d3DR,
# NewtonImpactRollingFrictionNSL, MoreauJeanOSI, RollingFrictionContact):
# one disk (R = 1, m = 1, I = 1/2, g = 9.81) on a line, the slope given by
# tilting the gravity. The vnative run puts the disk on a tilted segment and
# its trajectory is read in the frame of the segment. The cases:
#
#   roll      horizontal, rolling without slip, mu_r > 0: decelerates, stops
#   coulomb   horizontal, rolling without slip, mu_r = 0: rolls forever
#   slip      horizontal, sliding without rotation: slips, then rolls
#   incline   at rest on a slope, tan alpha > mu_r / R: rolls down
#   rest      at rest on a slope, tan alpha < mu_r / R: stays at rest
#
# The hdf5 tables are float32: the tolerance is on that scale. Each vnative
# run is in its own process (a second simulation in the same process
# crashes the vnative backend).
#
# usage: python3 rolling_friction_aos.py <disk_line_rolling_aos>
# exit status 0 if every case agrees

import math
import os
import subprocess
import sys
import tempfile

import h5py
import numpy as np

mu = 0.5
h = 1e-3
R, m, g = 1.0, 1.0, 9.81

# the cases: (mu_r, tan alpha, v0, w0, T)
CASES = {
    "roll": (0.1, 0.0, 2.0, -2.0, 4.0),
    "coulomb": (0.0, 0.0, 2.0, -2.0, 1.0),
    "slip": (0.1, 0.0, 2.0, 0.0, 2.0),
    "incline": (0.1, 0.2, 0.0, 0.0, 1.0),
    "rest": (0.1, 0.05, 0.0, 0.0, 0.5),
}


def run_vnative(case, filename):
    import siconos.numerics as sn
    from siconos.io.mechanics_run import (
        MechanicsHdf5Runner,
        MechanicsHdf5Runner_run_options,
        RunnerConfig,
    )
    from siconos.mechanics.collision.tools import Contactor

    import nonos  # noqa: F401 registers the vnative backend

    mu_r, tan, v0, w0, T = CASES[case]
    alpha = math.atan(tan)
    c, s = math.cos(alpha), math.sin(alpha)
    # the reference frame tilted by -alpha: the slope goes down along +x
    L = 50.0
    # exactly tangent to a tilted segment, the round-off of the gap can leave
    # the contact inactive for one step (a free fall of g h): the disk starts
    # 1e-9 inside the segment, below the float32 resolution of the tables
    d = R - 1e-9 if tan else R
    config = RunnerConfig("vnative")
    with MechanicsHdf5Runner(config=config, io_filename=filename) as io:
        io.add_primitive_shape("Disk", "Disk", [R])
        io.add_primitive_shape("Ground", "Segment", (-L * c, L * s, L * c, -L * s))
        io.add_Newton_impact_rolling_friction_nsl("contact", mu=mu, mu_r=mu_r, e=0.0)
        io.add_object(
            "disk",
            [Contactor("Disk")],
            translation=[d * s, d * c],
            orientation=[0],
            velocity=[v0 * c, -v0 * s, w0],
            mass=m,
            inertia=0.5 * m * R**2,
        )
        io.add_object("ground", [Contactor("Ground")], translation=[0, 0])

    options = sn.SolverOptions(sn.solver_ids.SICONOS_ROLLING_FRICTION_2D_NSGS)
    options.iparam[sn.params.SICONOS_IPARAM_MAX_ITER] = 1000
    options.dparam[sn.params.SICONOS_DPARAM_TOL] = 1e-12
    run_options = MechanicsHdf5Runner_run_options()
    run_options["t0"] = 0
    run_options["T"] = T
    run_options["h"] = h
    run_options["theta"] = 0.5
    run_options["solver_options"] = options
    # the runner divides its gravity, 9.80665, by gravity_scale
    run_options["gravity_scale"] = 9.80665 / g
    with MechanicsHdf5Runner(mode="r+", config=config, io_filename=filename) as io:
        io.run(run_options)


def vnative(case, tmpdir):
    """t, then x, vx and w in the frame of the line"""
    filename = os.path.join(tmpdir, case + ".hdf5")
    r = subprocess.run(
        [sys.executable, __file__, "--vnative", case, filename],
        capture_output=True,
        text=True,
    )
    if r.returncode != 0:
        raise RuntimeError(f"exit status {r.returncode}\n{r.stdout}{r.stderr}")
    with h5py.File(filename, "r") as f:
        q, v = f["data/dynamic"][:], f["data/velocities"][:]
    ds = v[0, 1]
    q, v = q[q[:, 1] == ds], v[v[:, 1] == ds]
    alpha = math.atan(CASES[case][1])
    c, s = math.cos(alpha), math.sin(alpha)
    return q[:, 0], q[:, 2] * c - q[:, 3] * s, v[:, 2] * c - v[:, 3] * s, v[:, 7]


def aos(program, case):
    """the columns t x y theta vx vy w rn rt rr"""
    mu_r, tan, v0, w0, T = CASES[case]
    args = [mu, mu_r, math.atan(tan), v0, w0, T, h, 0.5]
    r = subprocess.run(
        [program] + [repr(float(a)) for a in args],
        capture_output=True,
        text=True,
        check=True,
    )
    return np.loadtxt(r.stdout.splitlines())


if len(sys.argv) == 4 and sys.argv[1] == "--vnative":
    run_vnative(sys.argv[2], sys.argv[3])
    sys.exit(0)

program = os.path.abspath(sys.argv[1])
tmpdir = tempfile.mkdtemp()
failed = False
for case in CASES:
    ref = aos(program, case)
    try:
        t, x, vx, w = vnative(case, tmpdir)
    except RuntimeError as e:
        print(f"FAILED: {case}: {e}")
        failed = True
        continue
    rows = np.rint(t / h).astype(int)
    if rows.max() >= len(ref):
        print(f"FAILED: {case}: {len(t)} vnative rows, {len(ref)} reference rows")
        failed = True
        continue
    r = ref[rows]
    ex = np.max(np.abs(x - r[:, 1]))
    ev = np.max(np.abs(vx - r[:, 4]))
    ew = np.max(np.abs(w - r[:, 6]))
    # float32 tables: a few ulps of the largest value
    tol_x = 1e-6 * max(1.0, np.max(np.abs(r[:, 1])))
    tol_v = 1e-6 * max(1.0, np.max(np.abs(r[:, 4])), np.max(np.abs(r[:, 6])))
    ok = ex <= tol_x and ev <= tol_v and ew <= tol_v
    print(
        f"{case}: {len(t)} rows, max |x - x_ref| {ex:.2e}, max |vx - vx_ref| {ev:.2e}, "
        f"max |w - w_ref| {ew:.2e}; x(T) {x[-1]:.7f} ref {r[-1, 1]:.7f}"
        f"{'' if ok else ' FAILED'}"
    )
    failed |= not ok

sys.exit(1 if failed else 0)
