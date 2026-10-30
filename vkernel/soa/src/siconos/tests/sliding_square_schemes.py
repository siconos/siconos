# Check: the sliding square with the Moreau-Jean and Anitescu-Tasora schemes
# of the vnative backend, through the runner.
#
# A unit square launched at v0 on a ground segment with friction mu: with
# Moreau-Jean it stays on the ground (y = 0.5) and stops at v0^2 / (2 mu g);
# with Anitescu-Tasora (nonos.bridge.AnitescuTasoraOSI, the cone
# complementarity) it slides with u_N = mu |u_T| and lifts off. The C++
# ctest vkernel_nscd_cases checks the numbers; this script checks the python
# path (the scheme chosen by run_options["osi"]).
#
# Each scheme runs in its own process: a second simulation in the same
# process crashes the vnative backend (a segmentation fault, whatever the
# schemes).
#
# usage: python3 sliding_square_schemes.py
# exit status 0 if Moreau-Jean stays on the ground and Anitescu-Tasora lifts off

import subprocess
import sys
import tempfile
import os

import h5py
import numpy as np

import siconos.numerics as sn

from siconos.mechanics.collision.tools import Contactor

from siconos.io.mechanics_run import (
    MechanicsHdf5Runner,
    MechanicsHdf5Runner_run_options,
    RunnerConfig,
)
import nonos  # noqa: F401 registers the vnative backend
import nonos.bridge
from siconos import integrators

v0 = 2.0
mu = 0.3
g = 9.81
h = 0.005


def run(osi, filename):
    runner_config = RunnerConfig("vnative")
    with MechanicsHdf5Runner(config=runner_config, io_filename=filename) as io:
        io.add_primitive_shape("Ground", "Segment", (-20, 0, 20, 0))
        io.add_object("ground", [Contactor("Ground")], translation=[0, 0])
        square = np.array([[-0.5, -0.5], [0.5, -0.5], [0.5, 0.5], [-0.5, 0.5]])
        io.add_convex_shape("Square", square)
        io.add_object(
            "square",
            [Contactor("Square")],
            translation=[0, 0.5],
            orientation=[0],
            velocity=[v0, 0, 0],
            mass=1,
            inertia=1 / 6,
        )
        io.add_Newton_impact_friction_nsl("contact", mu=mu, e=0)

    options = sn.SolverOptions(sn.solver_ids.SICONOS_FRICTION_2D_NSGS)
    options.iparam[sn.params.SICONOS_IPARAM_MAX_ITER] = 100000
    options.dparam[sn.params.SICONOS_DPARAM_TOL] = 1e-12

    run_options = MechanicsHdf5Runner_run_options()
    run_options["t0"] = 0
    run_options["T"] = 0.5
    run_options["h"] = h
    run_options["osi"] = osi
    run_options["solver_options"] = options
    run_options["output_frequency"] = None
    with MechanicsHdf5Runner(
        mode="r+", config=runner_config, io_filename=filename
    ) as io:
        io.run(run_options)

    # rows (t, id, x, y, theta, ...) of data/dynamic: the square is the only
    # dynamic object
    with h5py.File(filename, "r") as f:
        dynamic = f["data/dynamic"][:]
    return dynamic[:, 3]


if len(sys.argv) == 3:
    # one scheme: the heights of the square, in a file
    osi = {"mj": integrators.MoreauJeanOSI, "at": nonos.bridge.AnitescuTasoraOSI}
    np.save(sys.argv[2], run(osi[sys.argv[1]], sys.argv[2] + ".hdf5"))
    sys.exit(0)

tmpdir = tempfile.mkdtemp()
heights = {}
for scheme in ("mj", "at"):
    out = os.path.join(tmpdir, scheme + ".npy")
    subprocess.run(
        [sys.executable, __file__, scheme, out], check=True, stdout=subprocess.DEVNULL
    )
    heights[scheme] = np.load(out)

mj_lift = np.abs(heights["mj"] - 0.5).max()
at_lift = (heights["at"] - 0.5).max()
print(f"Moreau-Jean: max |y - 0.5| = {mj_lift:.3e}")
print(f"Anitescu-Tasora: max lift-off = {at_lift:.3e} (h mu v0 = {h * mu * v0:.3e})")

# Moreau-Jean: no lift-off, but the vnative backend keeps a constraint
# activation threshold of 0 (run_options["constraint_activation_threshold"]
# goes to the kernel integrator object, not to vnative), so a contact resting
# at gap 0 may drop out of the active set on a rounding error and the square
# falls by ~ g h^2 on its other corner (vkernel_nscd_cases uses 1e-10)
if mj_lift > 2 * g * h * h:
    print("FAILED: the Moreau-Jean square left the ground")
    sys.exit(1)
if at_lift < 0.1 * h * mu * v0:
    print("FAILED: the Anitescu-Tasora square did not lift off")
    sys.exit(1)
print("OK")
