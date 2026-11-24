# Siconos is a program dedicated to modeling, simulation and control
# of non smooth dynamical systems.
#
# SPDX-FileCopyrightText: INRIA
# SPDX-License-Identifier: Apache-2.0
#
# Check: the quasi inelastic shock law of Jean 1999 (Eq 25) on the vnative
# backend, through the runner.
#
# A disk (R = 0.1, m = 1, I = 1/2 m R^2) is dropped on a horizontal floor,
# v = (0, -1), e = 0, mu = 1.  The law replaces the restitution: it closes the
# gap over the step and never ejects.  The cases:
#
#   qi-half   MoreauJeanQuasiInelasticOSI, theta = 0.5   no ejection, rest
#   qi-one    MoreauJeanQuasiInelasticOSI, theta = 1     gap closes exactly
#   newton    MoreauJeanOSI, theta = 0.5, e = 0          the comparison
#   e-positive  the quasi inelastic law with e > 0        an error
#   fremond   the quasi inelastic law and the Fremond law  an error
#
# The two-step property: the cumulated normal impulse of the cf table equals
# the single-step e = 0 impulse m |v_n^-| (Jean sec 4.2, positive homogeneity),
# and it is spread over more than one step (the gap is closed first).
#
# Each case runs in its own process: a second simulation in the same process
# crashes the vnative backend.
#
# usage: python3 quasi_inelastic.py
# exit status 0 if every case gives its result

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
v0 = [0.0, -1.0, 0.0]
y0 = R + 0.02


def run(case, filename):
    config = RunnerConfig("vnative")
    with MechanicsHdf5Runner(config=config, io_filename=filename) as io:
        io.add_primitive_shape("Disk", "Disk", [R])
        if case.startswith("qi-incline"):
            # a disk on an incline (gravity on): a disk ROLLS (friction prevents
            # slipping, not rolling).  mu high: rolls without slip; mu low: slips.
            import math

            alpha = math.radians(20.0)
            mu_law = 0.6 if case == "qi-incline-roll" else 0.02
            L = 4.0
            dx, dy = math.cos(alpha), math.sin(alpha)
            nx, ny = -math.sin(alpha), math.cos(alpha)  # normal, material side down
            io.add_primitive_shape("Ground", "Segment", (0, 0, L * dx, L * dy))
            cx, cy = 0.5 * L * dx + R * nx, 0.5 * L * dy + R * ny
            io.add_object(
                "disk",
                [Contactor("Disk")],
                translation=[cx, cy],
                orientation=[0],
                velocity=[0.0, 0.0, 0.0],
                mass=m,
                inertia=inertia,
            )
            io.add_object("ground", [Contactor("Ground")], translation=[0, 0])
            io.add_Newton_impact_friction_nsl("contact", mu=mu_law, e=0.0)
        else:
            io.add_primitive_shape("Ground", "Segment", (-10, 0, 10, 0))
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
            e = 1.0 if case == "e-positive" else 0.0
            if case == "fremond":
                io.add_Fremond_impact_friction_nsl("contact", mu=1.0, e=0.0)
            else:
                io.add_Newton_impact_friction_nsl("contact", mu=1.0, e=e)

    options = sn.SolverOptions(sn.solver_ids.SICONOS_FRICTION_2D_NSGS)
    options.iparam[sn.params.SICONOS_IPARAM_MAX_ITER] = 1000
    options.dparam[sn.params.SICONOS_DPARAM_TOL] = 1e-12

    run_options = MechanicsHdf5Runner_run_options()
    run_options["t0"] = 0.0
    run_options["T"] = 0.5 if case.startswith("qi-incline") else 0.3
    run_options["h"] = 1e-3
    run_options["theta"] = 1.0 if case == "qi-one" else 0.5
    run_options["solver_options"] = options
    if case in ("qi-half", "qi-one", "e-positive", "fremond") or case.startswith(
        "qi-incline"
    ):
        run_options["osi"] = nonos.bridge.MoreauJeanQuasiInelasticOSI
    with MechanicsHdf5Runner(mode="r+", config=config, io_filename=filename) as io:
        io.run(run_options)

    with h5py.File(filename, "r") as f:
        dyn = f["data/dynamic"][:]
        vel = f["data/velocities"][:]
        cf = f["data/cf"][:]
    return dyn, vel, cf


if len(sys.argv) == 3:
    case, out = sys.argv[1:]
    try:
        dyn, vel, cf = run(case, out + ".hdf5")
        np.save(out + ".dyn.npy", dyn)
        np.save(out + ".vel.npy", vel)
        np.save(out + ".cf.npy", cf)
    except ValueError as e:
        print(f"ValueError: {e}")
        sys.exit(3)
    sys.exit(0)


def load(out, suffix):
    return np.load(out + suffix)


failed = False
tmpdir = tempfile.mkdtemp()
ptol = 1e-4  # an upward normal velocity allowed after the landing
summary = {}


def analyse(case):
    """the disk trajectory and the cumulated normal impulse of a run"""
    out = os.path.join(tmpdir, case + ".npy")
    r = subprocess.run(
        [sys.executable, __file__, case, out], capture_output=True, text=True
    )
    if r.returncode == 3:
        return None, r.stdout.strip()
    if r.returncode != 0:
        print(f"FAILED: {case} exited {r.returncode}\n{r.stderr}")
        sys.exit(1)
    dyn = load(out, ".dyn.npy")
    vel = load(out, ".vel.npy")
    cf = load(out, ".cf.npy")
    did = [i for i in np.unique(dyn[:, 1]) if i > 0][0]
    d = dyn[dyn[:, 1] == did]
    v = vel[vel[:, 1] == did]
    return (d[:, 3], v[:, 3], cf), None


for case in ("qi-half", "qi-one", "newton"):
    res, err = analyse(case)
    y, vy, cf = res
    landed = int(np.argmax(y <= R + 2e-3))
    vy_after = vy[landed:]
    vb = vy[max(landed - 1, 0)]
    # the impact impulse: the cumulated normal impulse r_n over the first three
    # output steps where the contact is active (the shock window), not the whole
    # resting phase (whose r_n supports gravity every step)
    if cf.shape[0]:
        times = np.unique(cf[:, 0])
        window = times[:4]
        pn = cf[np.isin(cf[:, 0], window), 20].sum()
        n_contact_steps = int(
            np.count_nonzero([cf[cf[:, 0] == t, 20].sum() > 1e-6 for t in window])
        )
    else:
        pn, n_contact_steps = 0.0, 0
    summary[case] = (vy_after.max(), y[-1] - R, vy[-1], pn, n_contact_steps, vb)
    print(
        f"{case}: max vy after landing={vy_after.max():.3e} "
        f"final y-R={y[-1] - R:.3e} final vy={vy[-1]:.3e} "
        f"Pn={pn:.4f} steps={n_contact_steps} v_before={vb:.4f}"
    )
    if case in ("qi-half", "qi-one"):
        if vy_after.max() > ptol:
            print(f"FAILED: {case} ejects (max vy after landing {vy_after.max():.3e})")
            failed = True
        if abs(vy[-1]) > 1e-3:
            print(f"FAILED: {case} does not come to rest (vy {vy[-1]:.3e})")
            failed = True
        # the two-step impulse: cumulated normal impulse over the shock window
        # == m |v_n^-| (Jean sec 4.2, positive homogeneity)
        if abs(pn - m * abs(vb)) > 0.10 * m * abs(vb):
            print(f"FAILED: {case} impulse {pn:.4f} != m|v-| {m * abs(vb):.4f}")
            failed = True
        if n_contact_steps < 2:
            print(f"FAILED: {case} the shock is not split over two steps")
            failed = True

# theta = 1 closes the gap exactly, theta = 0.5 leaves the documented
# permanent overlap of g_pos (no push-out)
if abs(summary["qi-one"][1]) > 1e-6:
    print(f"FAILED: qi-one gap not closed ({summary['qi-one'][1]:.3e})")
    failed = True
if not (-2e-3 < summary["qi-half"][1] < 1e-6):
    print(f"FAILED: qi-half overlap out of range ({summary['qi-half'][1]:.3e})")
    failed = True

# the incline (gravity on): a disk ROLLS (friction stops the slipping, not the
# rolling).  mu high: rolls without slip; mu low: slips.  Measured on the
# contact tangential slip  v_t - omega R  along the incline.
import math

ALPHA = math.radians(20.0)
for case in ("qi-incline-roll", "qi-incline-slip"):
    out = os.path.join(tmpdir, case + ".npy")
    r = subprocess.run(
        [sys.executable, __file__, case, out], capture_output=True, text=True
    )
    if r.returncode != 0:
        print(f"FAILED: {case} exited {r.returncode}\n{r.stderr[-800:]}")
        failed = True
        continue
    dyn = load(out, ".dyn.npy")
    vel = load(out, ".vel.npy")
    did = [i for i in np.unique(dyn[:, 1]) if i > 0][0]
    d = dyn[dyn[:, 1] == did]
    v = vel[vel[:, 1] == did]
    disp = float(np.hypot(d[-1, 2] - d[0, 2], d[-1, 3] - d[0, 3]))
    vt = v[-1, 2] * math.cos(ALPHA) + v[-1, 3] * math.sin(ALPHA)
    slip = float(abs(vt + v[-1, 7] * R))
    print(f"{case}: displacement={disp:.3e} tangential slip={slip:.3e}")
    if disp < 2e-3:
        print(f"FAILED: {case} does not move ({disp:.3e})")
        failed = True
    if case == "qi-incline-roll" and slip > 1e-2:
        print(f"FAILED: {case} slips (slip {slip:.3e})")
        failed = True
    if case == "qi-incline-slip" and slip < 0.05:
        print(f"FAILED: {case} does not slip (slip {slip:.3e})")
        failed = True

# the refusals
for case, what in (("e-positive", "e = 0"), ("fremond", "Frémond")):
    out = os.path.join(tmpdir, case + ".npy")
    r = subprocess.run(
        [sys.executable, __file__, case, out], capture_output=True, text=True
    )
    print(f"{case}: exit {r.returncode}, {r.stdout.strip()}")
    if r.returncode != 3:
        print(f"FAILED: {case} should be refused")
        failed = True
    elif case == "e-positive" and "replaces the restitution" not in r.stdout:
        print(f"FAILED: {case} wrong message")
        failed = True
    elif case == "fremond" and "cannot be mixed" not in r.stdout:
        print(f"FAILED: {case} wrong message")
        failed = True

sys.exit(1 if failed else 0)
