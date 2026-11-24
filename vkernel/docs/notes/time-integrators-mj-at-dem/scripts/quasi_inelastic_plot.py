# Siconos is a program dedicated to modeling, simulation and control
# of non smooth dynamical systems.
#
# SPDX-FileCopyrightText: INRIA
# SPDX-License-Identifier: Apache-2.0
#
# The comparison figure of the Jean quasi inelastic shock law, on the same
# drop and the same incline (R = 0.1, m = 1, e = 0, theta = 0.5):
#   - the Newton impact law (plain Moreau-Jean, velocity Signorini),
#   - Jean 1999 Eq 25 (MoreauJeanQuasiInelasticOSI, gap -> 0, no push-out),
#   - the Anitescu-Tasora stabilisation (AnitescuTasoraInelastic, signed gap
#     clamped at the recovery speed: an overlap is pushed out).
#
# Top row, a disk dropped on a floor: the gap y - R and the normal velocity.
# The quasi inelastic law closes the gap without ejecting, the AT stabilisation
# pushes the overlap out, Newton stops on the first active step.
# Bottom row, a disk on a 20 deg incline, mu = 0.6: it rolls (friction stops the
# slip, not the rolling) - the gap stays at 0 under all three laws.
#
# usage: python3 quasi_inelastic_plot.py [out.png]
# Each case runs in its own process (a second simulation in the same process
# crashes the vnative backend).

import math
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
ALPHA = math.radians(20.0)

LAWS = {
    "newton": ("MoreauJean", "Newton, e = 0", dict(color="0.55", ls="--", lw=1.4)),
    "quasi-inelastic": (
        "MoreauJeanQuasiInelastic",
        "Jean 1999 Eq 25",
        dict(color="C0", ls="-", lw=1.8),
    ),
    "at-inelastic": (
        "AnitescuTasoraInelastic",
        "Anitescu-Tasora stabilisation",
        dict(color="C3", ls="-.", lw=1.4),
    ),
}


def run(case, law, scene, filename):
    config = RunnerConfig("vnative")
    with MechanicsHdf5Runner(config=config, io_filename=filename) as io:
        io.add_primitive_shape("Disk", "Disk", [R])
        if scene == "drop":
            io.add_primitive_shape("Ground", "Segment", (-10, 0, 10, 0))
            io.add_object(
                "disk",
                [Contactor("Disk")],
                translation=[0, R + 0.02],
                orientation=[0],
                velocity=[0.0, -1.0, 0.0],
                mass=m,
                inertia=inertia,
            )
            io.add_object("ground", [Contactor("Ground")], translation=[0, 0])
            io.add_Newton_impact_friction_nsl("contact", mu=1.0, e=0.0)
        else:  # incline
            dx, dy = math.cos(ALPHA), math.sin(ALPHA)
            nx, ny = -math.sin(ALPHA), math.cos(ALPHA)
            io.add_primitive_shape("Ground", "Segment", (0, 0, 4 * dx, 4 * dy))
            cx, cy = 2.0 * dx + R * nx, 2.0 * dy + R * ny
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
            io.add_Newton_impact_friction_nsl("contact", mu=0.6, e=0.0)

    options = sn.SolverOptions(sn.solver_ids.SICONOS_FRICTION_2D_NSGS)
    options.iparam[sn.params.SICONOS_IPARAM_MAX_ITER] = 1000
    options.dparam[sn.params.SICONOS_DPARAM_TOL] = 1e-10

    ro = MechanicsHdf5Runner_run_options()
    ro["t0"] = 0.0
    ro["T"] = 0.04 if scene == "drop" else 0.5
    ro["h"] = 1e-3
    ro["theta"] = 0.5
    ro["solver_options"] = options
    if law in ("MoreauJeanQuasiInelastic", "AnitescuTasoraInelastic"):
        ro["osi"] = getattr(
            nonos.bridge,
            law.replace("MoreauJeanQuasiInelastic", "MoreauJeanQuasiInelasticOSI"),
        )
    with MechanicsHdf5Runner(mode="r+", config=config, io_filename=filename) as io:
        io.run(ro)

    with h5py.File(filename, "r") as f:
        dyn = f["data/dynamic"][:]
        vel = f["data/velocities"][:]
    did = [i for i in np.unique(dyn[:, 1]) if i > 0][0]
    d = dyn[dyn[:, 1] == did]
    v = vel[vel[:, 1] == did]
    if scene == "drop":
        return d[:, 0], d[:, 3] - R, v[:, 3], None
    vt = v[:, 2] * math.cos(ALPHA) + v[:, 3] * math.sin(ALPHA)
    slip = np.abs(vt + v[:, 7] * R)
    # the signed distance to the incline line
    g = (d[:, 2] * math.sin(ALPHA) - d[:, 3] * math.cos(ALPHA)) + R
    dist = (d[:, 2] - d[0, 2]) * math.cos(ALPHA) + (d[:, 3] - d[0, 3]) * math.sin(ALPHA)
    return d[:, 0], g, dist, slip


if len(sys.argv) == 5:
    case, law, scene, out = sys.argv[1:]
    t, a, b, c = run(case, law, scene, out + ".hdf5")
    np.savez(out + ".npz", t=t, a=a, b=b, c=c)
    sys.exit(0)


out_png = sys.argv[1] if len(sys.argv) > 1 else "quasi_inelastic.png"
tmpdir = tempfile.mkdtemp()
data = {}
for scene in ("drop", "incline"):
    for name, (law, _, _) in LAWS.items():
        out = os.path.join(tmpdir, f"{scene}_{name}")
        subprocess.run([sys.executable, __file__, name, law, scene, out], check=True)
        data[(scene, name)] = np.load(out + ".npz")

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402

fig, axes = plt.subplots(2, 2, figsize=(9.0, 6.2))
for name, (_, label, style) in LAWS.items():
    d = data[("drop", name)]
    axes[0, 0].plot(d["t"], d["a"], label=label, **style)
    axes[0, 1].plot(d["t"], d["b"], label=label, **style)
    di = data[("incline", name)]
    axes[1, 0].plot(di["t"], di["a"], label=label, **style)
    axes[1, 1].plot(di["t"], di["b"], label=label, **style)

axes[0, 0].set(ylabel="gap y-R [m]", title="drop: gap (overlap < 0)")
axes[0, 1].set(ylabel="normal velocity $v_y$ [m/s]", title="drop: normal velocity")
axes[1, 0].set(xlabel="t [s]", ylabel="gap to incline [m]", title="incline: gap")
axes[1, 1].set(
    xlabel="t [s]",
    ylabel="distance along slope [m]",
    title="incline: rolling distance",
)
axes[0, 0].set_xlim(0.012, 0.030)
axes[1, 0].set_ylim(-5e-4, 5e-4)
for ax in axes.ravel():
    ax.axhline(0, color="k", lw=0.5)
    ax.grid(alpha=0.3)
axes[0, 0].legend(fontsize=7, loc="lower left")
fig.tight_layout()
fig.savefig(out_png, dpi=130)
print(f"wrote {out_png}")
