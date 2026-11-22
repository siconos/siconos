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
# The comparison figure of the Jean quasi inelastic shock law: a disk (R = 0.1,
# m = 1) dropped on a floor at v = (0, -1), e = 0, mu = 1, theta = 0.5, under
#   - the Newton impact law (plain Moreau-Jean, velocity Signorini),
#   - Jean 1999 Eq 25 (MoreauJeanQuasiInelasticOSI, gap -> 0, no push-out),
#   - the Anitescu-Tasora stabilisation (AnitescuTasoraInelastic, signed gap
#     clamped at the recovery speed: an overlap is pushed out).
#
# The gap y - R and the normal velocity vy around the landing show the two
# mechanisms: the quasi inelastic law closes the gap and never ejects, the AT
# stabilisation pushes the overlap out, the Newton law stops on the first
# active step.
#
# usage: python3 quasi_inelastic_plot.py [out.png]
# Each case runs in its own process (a second simulation in the same process
# crashes the vnative backend).

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
y0 = R + 0.02

CASES = {
    "newton": ("MoreauJean", "Newton, e = 0"),
    "quasi-inelastic": ("MoreauJeanQuasiInelastic", "Jean 1999 Eq 25"),
    "at-inelastic": ("AnitescuTasoraInelastic", "Anitescu-Tasora stabilisation"),
}


def run(case, filename):
    config = RunnerConfig("vnative")
    with MechanicsHdf5Runner(config=config, io_filename=filename) as io:
        io.add_primitive_shape("Disk", "Disk", [R])
        io.add_primitive_shape("Ground", "Segment", (-10, 0, 10, 0))
        io.add_object(
            "disk",
            [Contactor("Disk")],
            translation=[0, y0],
            orientation=[0],
            velocity=[0.0, -1.0, 0.0],
            mass=m,
            inertia=inertia,
        )
        io.add_object("ground", [Contactor("Ground")], translation=[0, 0])
        io.add_Newton_impact_friction_nsl("contact", mu=1.0, e=0.0)

    options = sn.SolverOptions(sn.solver_ids.SICONOS_FRICTION_2D_NSGS)
    options.iparam[sn.params.SICONOS_IPARAM_MAX_ITER] = 1000
    options.dparam[sn.params.SICONOS_DPARAM_TOL] = 1e-12

    ro = MechanicsHdf5Runner_run_options()
    ro["t0"] = 0.0
    ro["T"] = 0.04
    ro["h"] = 1e-3
    ro["theta"] = 0.5
    ro["solver_options"] = options
    if case == "quasi-inelastic":
        ro["osi"] = nonos.bridge.MoreauJeanQuasiInelasticOSI
    elif case == "at-inelastic":
        ro["osi"] = nonos.bridge.AnitescuTasoraInelastic
    with MechanicsHdf5Runner(mode="r+", config=config, io_filename=filename) as io:
        io.run(ro)

    with h5py.File(filename, "r") as f:
        dyn = f["data/dynamic"][:]
        vel = f["data/velocities"][:]
    did = [i for i in np.unique(dyn[:, 1]) if i > 0][0]
    d = dyn[dyn[:, 1] == did]
    v = vel[vel[:, 1] == did]
    return d[:, 0], d[:, 3], v[:, 3]


if len(sys.argv) == 3:
    case, out = sys.argv[1:]
    t, y, vy = run(case, out + ".hdf5")
    np.savez(out + ".npz", t=t, y=y, vy=vy)
    sys.exit(0)


out_png = sys.argv[1] if len(sys.argv) > 1 else "quasi_inelastic.png"
tmpdir = tempfile.mkdtemp()
data = {}
for case in CASES:
    out = os.path.join(tmpdir, case)
    subprocess.run([sys.executable, __file__, case, out], check=True)
    data[case] = np.load(out + ".npz")

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402

fig, (ax_y, ax_v) = plt.subplots(1, 2, figsize=(9.0, 3.4))
styles = {
    "newton": dict(color="0.55", ls="--", lw=1.4),
    "quasi-inelastic": dict(color="C0", ls="-", lw=1.8),
    "at-inelastic": dict(color="C3", ls="-.", lw=1.4),
}
for case, (_, label) in CASES.items():
    d = data[case]
    ax_y.plot(d["t"], d["y"] - R, label=label, **styles[case])
    ax_v.plot(d["t"], d["vy"], label=label, **styles[case])

ax_y.set(xlabel="t [s]", ylabel="gap y - R [m]", title="gap (overlap < 0)")
ax_v.set(
    xlabel="t [s]",
    ylabel="normal velocity $v_y$ [m/s]",
    title="normal velocity (ejection > 0)",
)
ax_y.axhline(0, color="k", lw=0.5)
ax_v.axhline(0, color="k", lw=0.5)
for ax in (ax_y, ax_v):
    ax.set_xlim(0.012, 0.030)
    ax.grid(alpha=0.3)
ax_y.legend(fontsize=7, loc="lower left")
fig.tight_layout()
fig.savefig(out_png, dpi=130)
print(f"wrote {out_png}")
