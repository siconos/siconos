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

# Check: rolling resistance (NewtonImpactRollingFrictionNSL) on the vnative
# backend, through the runner, against the analytical values.
#
# A disk (R, m, I = m R^2 / 2) rolls without slip on a segment, the rolling
# resistance mu_r a length (Acary & Bourrier 2021): |m_r| <= mu_r r_n. On a
# slope alpha, rolling without slip:
#
#   a = (2/3) g (sin alpha - mu_r / R cos alpha)
#
# so that on the horizontal (v0 > 0) it decelerates by (2/3) g mu_r / R and
# stops after 3 v0^2 R / (4 g mu_r); with mu_r = 0 it rolls forever; at rest
# on a slope it stays at rest iff tan alpha <= mu_r / R. The cases:
#
#   roll           horizontal, mu_r > 0, Moreau-Jean: deceleration, stop
#   roll-mjf       the same, the Moreau-Jean-Frémond scheme
#   coulomb        horizontal, mu_r = 0: rolls forever
#   incline        at rest on slopes of tan alpha (scan around mu_r / R),
#                  one disk per slope, side by side in one scene
#   roll-at        horizontal, Anitescu-Tasora (projected Gauss-Seidel): the
#                  associated (convexified) law, the disk decelerates
#   roll-apgd      the same, APGD
#   pair           two disks rolling on the ground against each other
#   pile           a small pile of disks on the ground
#   gosi, dem      rolling resistance with MoreauJeanGOSI, ExplicitDEMOSI:
#                  refused
#
# The Moreau-Jean-Frémond scheme stops with an alternating velocity below
# a h. roll, coulomb, roll-at, roll-apgd, pair and pile are also checked against
# their reference trajectories (ref_trajectory.py). Each case runs in its own
# process: a second simulation in the same process crashes the vnative
# backend.
#
# usage: python3 rolling_friction.py
# exit status 0 if every case gives its value

import functools
import math
import os
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

# g: the gravity of the runner
R, m, g = 0.5, 1.0, 9.80665
inertia = 0.5 * m * R**2
mu, mu_r = 0.5, 0.05
v0 = 1.0
# the checks hold at any h (the deceleration is exact, the tolerances of the
# stop scale with h): a coarse one keeps the test short
h = 1e-2

# the cases: (T, mu_r, tan alpha)
CASES = {
    "roll": (2.0, mu_r, 0.0),
    "roll-mjf": (2.0, mu_r, 0.0),
    "coulomb": (1.0, 0.0, 0.0),
    "roll-at": (1.0, mu_r, 0.0),
    "roll-apgd": (1.0, mu_r, 0.0),
    "pair": (1.0, mu_r, 0.0),
    "pile": (1.0, mu_r, 0.0),
    "gosi": (0.1, mu_r, 0.0),
    "dem": (0.1, mu_r, 0.0),
}
# the slopes of the incline case, SPACING apart
TANS = (0.09, 0.099, 0.101, 0.11, 0.2)
SPACING = 20.0
CASES["incline"] = (0.5, mu_r, 0.0)


def run(case, filename):
    T, case_mu_r, tan = CASES[case]
    alpha = math.atan(tan)
    rolling = case.startswith(("roll", "coulomb", "pair"))
    config = RunnerConfig("vnative")
    with MechanicsHdf5Runner(config=config, io_filename=filename) as io:
        io.add_primitive_shape("Disk", "Disk", [R])
        # the ground through the origin, of slope tan alpha
        L = 50.0
        io.add_primitive_shape("Ground", "Segment", (-L, -L * tan, L, L * tan))
        n = (-math.sin(alpha), math.cos(alpha))
        if case == "incline":
            # a slope through (k SPACING, 0) and a disk at rest on it, each
            for k, slope in enumerate(TANS):
                x0, a = k * SPACING, math.atan(slope)
                io.add_primitive_shape(
                    f"Slope{k}", "Segment", (x0 - 5, -5 * slope, x0 + 5, 5 * slope)
                )
                io.add_object(f"slope{k}", [Contactor(f"Slope{k}")], translation=[0, 0])
                io.add_object(
                    f"disk{k}",
                    [Contactor("Disk")],
                    translation=[x0 - R * math.sin(a), R * math.cos(a)],
                    orientation=[0],
                    mass=m,
                    inertia=inertia,
                )
        elif case == "pair":
            # side by side, the first one rolling into the second one
            for k, (x, v) in enumerate([(0.0, v0), (2 * R + 0.1, 0.0)]):
                io.add_object(
                    f"disk{k}",
                    [Contactor("Disk")],
                    translation=[x, R],
                    orientation=[0],
                    velocity=[v, 0, -v / R],
                    mass=m,
                    inertia=inertia,
                )
        elif case == "pile":
            # three rows, a triangle
            k = 0
            for row, count in enumerate((4, 3, 2)):
                for i in range(count):
                    x = (2 * i + row) * R * 1.001
                    y = R + row * R * math.sqrt(3) * 1.001
                    io.add_object(
                        f"disk{k}",
                        [Contactor("Disk")],
                        translation=[x, y],
                        orientation=[0],
                        mass=m,
                        inertia=inertia,
                    )
                    k += 1
        else:
            # on the ground, rolling without slip
            v = v0 if rolling else 0.0
            io.add_object(
                "disk",
                [Contactor("Disk")],
                translation=[R * n[0], R * n[1]],
                orientation=[0],
                velocity=[v, 0, -v / R],
                mass=m,
                inertia=inertia,
            )
        if case != "incline":
            io.add_object("ground", [Contactor("Ground")], translation=[0, 0])
        io.add_Newton_impact_rolling_friction_nsl(
            "contact", mu=mu, mu_r=case_mu_r, e=0.0
        )

    if case == "roll-mjf":
        from nonos import simoptions

        simoptions.apply({"scheme": "MoreauJeanFremond"})

    # a rolling solver id, mapped to the 2D friction one by the bridge
    options = sn.SolverOptions(sn.solver_ids.SICONOS_ROLLING_FRICTION_2D_NSGS)
    options.iparam[sn.params.SICONOS_IPARAM_MAX_ITER] = 1000
    options.dparam[sn.params.SICONOS_DPARAM_TOL] = 1e-12

    run_options = MechanicsHdf5Runner_run_options()
    run_options["t0"] = 0
    run_options["T"] = T
    run_options["h"] = h
    run_options["theta"] = 0.5
    run_options["solver_options"] = options
    run_options["output_frequency"] = 10
    if case == "roll-at":
        run_options["osi"] = nonos.bridge.AnitescuTasoraOSI
    elif case == "roll-apgd":
        run_options["osi"] = nonos.bridge.AnitescuTasoraAPGDOSI
    elif case == "gosi":
        from siconos import integrators

        run_options["osi"] = integrators.MoreauJeanGOSI
    elif case == "dem":
        run_options["osi"] = functools.partial(
            nonos.bridge.ExplicitDEMOSI, normal_stiffness=1e5
        )
    with MechanicsHdf5Runner(mode="r+", config=config, io_filename=filename) as io:
        io.run(run_options)

    # rows (t, id, vx, vy, vz, wx, wy, wz) of data/velocities, (t, id, x, y,
    # ...) of data/dynamic
    with h5py.File(filename, "r") as f:
        return f["data/velocities"][:], f["data/dynamic"][:]


if len(sys.argv) == 3:
    case, out = sys.argv[1:]
    try:
        v, q = run(case, out + ".hdf5")
        np.savez(out, v=v, q=q)
    except ValueError as e:
        print(f"ValueError: {e}")
        sys.exit(3)
    sys.exit(0)

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import case_processes  # noqa: E402
import ref_trajectory  # noqa: E402

failed = False


def fail(message):
    global failed
    print(f"FAILED: {message}")
    failed = True


tmpdir = tempfile.mkdtemp()
# all the cases side by side, then checked in turn
outputs = {case: os.path.join(tmpdir, case + ".npz") for case in CASES}
results = dict(
    zip(
        CASES,
        case_processes.run(
            [[sys.executable, __file__, case, outputs[case]] for case in CASES]
        ),
    )
)


def run_case(case):
    return results[case], outputs[case]


def disk_rows(out):
    data = np.load(out)
    v, q = data["v"], data["q"]
    ds = v[0, 1]
    return v[v[:, 1] == ds], q[q[:, 1] == ds]


# the hdf5 tables are float32
tol = 1e-5
a_roll = 2 / 3 * g * mu_r / R
t_stop = v0 / a_roll
x_stop = 3 * v0**2 * R / (4 * g * mu_r)

for case in ("roll", "roll-mjf"):
    r, out = run_case(case)
    if r.returncode != 0:
        fail(f"{case}: exit status {r.returncode}\n{r.stdout}{r.stderr}")
        continue
    v, q = disk_rows(out)
    t, vx, w, x = v[:, 0], v[:, 2], v[:, 7], q[:, 2]
    rolling = t < t_stop - 0.05
    a = -np.polyfit(t[rolling & (t > 0.05)], vx[rolling & (t > 0.05)], 1)[0]
    slip = np.max(np.abs(vx + w * R))
    print(
        f"{case}: deceleration {a:.6f} (expected {a_roll:.6f}), stop at x = {x[-1]:.6f} "
        f"(expected {x_stop:.6f}), vx(T) = {vx[-1]:.2e}, slip {slip:.1e}"
    )
    if abs(a - a_roll) > 1e-4 * a_roll:
        fail(f"{case}: deceleration")
    # Moreau-Jean-Frémond: at rest, the friction on u(k + 1/2) leaves an
    # alternating velocity below a h
    if abs(x[-1] - x_stop) > 2 * v0 * h or abs(vx[-1]) > (
        tol if case == "roll" else a_roll * h
    ):
        fail(f"{case}: stop")
    if slip > tol:
        fail(f"{case}: slip")

r, out = run_case("coulomb")
if r.returncode != 0:
    fail(f"coulomb: exit status {r.returncode}\n{r.stdout}{r.stderr}")
else:
    v, q = disk_rows(out)
    print(
        f"coulomb: vx(T) = {v[-1, 2]:.9f} (expected {v0}), slip {abs(v[-1, 2] + v[-1, 7] * R):.1e}"
    )
    if np.max(np.abs(v[:, 2] - v0)) > tol:
        fail("coulomb: the disk does not roll forever")
    if not ref_trajectory.check(out + ".hdf5", "rolling_friction_coulomb", 1e-8):
        fail("coulomb: reference")

# at rest on a slope: rolls iff tan alpha > mu_r / R; the disk of a slope
# by its first position
r, out = run_case("incline")
if r.returncode != 0:
    fail(f"incline: exit status {r.returncode}\n{r.stdout}{r.stderr}")
else:
    data = np.load(out)
    v, q = data["v"], data["q"]
    first = q[q[:, 0] == q[0, 0]]
    T = CASES["incline"][0]
    for k, tan in enumerate(TANS):
        ds = first[np.argmin(np.abs(first[:, 2] - k * SPACING)), 1]
        vk = v[v[:, 1] == ds]
        alpha = math.atan(tan)
        a = max(0.0, 2 / 3 * g * (math.sin(alpha) - mu_r / R * math.cos(alpha)))
        speed = math.hypot(vk[-1, 2], vk[-1, 3])
        print(f"incline-{tan}: speed(T) = {speed:.6e} (expected {a * T:.6e})")
        if abs(speed - a * T) > tol * max(1, a * T):
            fail(f"incline-{tan}: speed")

# Anitescu-Tasora: the associated law, the disk decelerates; reference only
for case in ("roll-at", "roll-apgd"):
    r, out = run_case(case)
    if r.returncode != 0:
        fail(f"{case}: exit status {r.returncode}\n{r.stdout}{r.stderr}")
        continue
    v, q = disk_rows(out)
    print(f"{case}: vx(T) = {v[-1, 2]:.6f}")
    if not v[-1, 2] < v0 - 0.1 * a_roll:
        fail(f"{case}: no deceleration")
    if not ref_trajectory.check(out + ".hdf5", f"rolling_friction_{case}", 1e-6):
        fail(f"{case}: reference")

for case in ("roll", "pair", "pile"):
    r, out = (
        run_case(case) if case != "roll" else (None, os.path.join(tmpdir, "roll.npz"))
    )
    if r is not None and r.returncode != 0:
        fail(f"{case}: exit status {r.returncode}\n{r.stdout}{r.stderr}")
        continue
    if not ref_trajectory.check(out + ".hdf5", f"rolling_friction_{case}", 1e-6):
        fail(f"{case}: reference")

# the schemes without rolling resistance refuse it
for case in ("gosi", "dem"):
    r, out = run_case(case)
    last = (r.stdout.strip().splitlines() or [""])[-1]
    print(f"{case}: exit status {r.returncode}, {last}")
    if r.returncode != 3 or "rolling resistance" not in r.stdout:
        fail(f"{case}: rolling resistance should be refused")

sys.exit(1 if failed else 0)
