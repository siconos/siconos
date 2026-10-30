# Check: two FEM squares in contact in the same steps.
#
# The FEM square of mesh_N5.msh and a copy of it shifted by 2 in x, a disk
# falling on each at the same place. The two halves of the scene are the
# same, so the two disks must follow the same vertical trajectory. Both FEM
# systems are involved at the same steps: the second one gets its velocity
# correction at its own dofs in the assembled vector.
#
# usage: python3 two_fem_bodies.py [--gosi | --dem] <mesh_N5.msh> [hdf5_output]
# --gosi: Moreau-Jean on the global problem (MoreauJeanGOSI, GFC2D_NSN_AC)
# --dem: the explicit DEM scheme (ExplicitDEMOSI, penalty springs; the FEM
#        systems keep the implicit theta method, the contact forces explicit)
# exit status 0 if the two disks follow the same trajectory

import functools
import os
import sys
import tempfile

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

gosi = "--gosi" in sys.argv[1:]
dem = "--dem" in sys.argv[1:]
args = [a for a in sys.argv[1:] if a not in ("--gosi", "--dem")]
mesh_filename = args[0]
hdf5_output = args[1] if len(args) > 1 else "two_fem_bodies.hdf5"

shift = 2.0
disk_radius = 0.05


def shifted_mesh(filename, dx):
    """the gmsh 2.2 file with its nodes shifted by dx"""
    lines = open(filename).read().splitlines()
    begin = lines.index("$Nodes")
    n = int(lines[begin + 1])
    for k in range(begin + 2, begin + 2 + n):
        i, x, y, z = lines[k].split()
        lines[k] = f"{i} {float(x) + dx:.16e} {y} {z}"
    return "\n".join(lines) + "\n"


tmpdir = tempfile.mkdtemp()
shifted_filename = os.path.join(tmpdir, "shifted.msh")
with open(shifted_filename, "w") as f:
    f.write(shifted_mesh(mesh_filename, shift))

positions = np.array(
    [[0.5, 1.0 + disk_radius + 1e-3], [0.5 + shift, 1.0 + disk_radius + 1e-3]]
)

runner_config = RunnerConfig("vnative")

with MechanicsHdf5Runner(config=runner_config, io_filename=hdf5_output) as io:
    io.add_primitive_shape("DiskR", "Disk", [disk_radius])
    velocities = np.zeros((2, 3), dtype=np.float64)
    velocities[:, 1] = -1.0
    io.add_objects(
        "disks",
        [Contactor("DiskR")],
        translations=positions,
        orientations=np.zeros(2, dtype=np.float64),
        velocities=velocities,
        mass=1,
        inertia=1,
    )

    io.add_primitive_shape("Ground-1", "Segment", (-10, 0, 10, 0))

    io.add_shape_data_from_file("Square", mesh_filename)
    io.add_object("square", [Contactor("Square")], material=(2500, 1e11, 0.3))
    io.add_shape_data_from_file("ShiftedSquare", shifted_filename)
    io.add_object(
        "shifted_square", [Contactor("ShiftedSquare")], material=(2500, 1e11, 0.3)
    )

    io.add_Newton_impact_friction_nsl("contact", mu=0.5, e=0)

if gosi:
    options = sn.SolverOptions(sn.solver_ids.SICONOS_GLOBAL_FRICTION_2D_NSN_AC)
else:
    options = sn.SolverOptions(sn.solver_ids.SICONOS_FRICTION_2D_NSGS)
options.iparam[sn.params.SICONOS_IPARAM_MAX_ITER] = 1000
options.dparam[sn.params.SICONOS_DPARAM_TOL] = 1e-12

run_options = MechanicsHdf5Runner_run_options()
run_options["t0"] = 0
run_options["T"] = 0.1
run_options["h"] = 0.005
run_options["theta"] = 0.50001
if gosi:
    run_options["osi"] = integrators.MoreauJeanGOSI
if dem:
    # h = 0.005 < 2 sqrt(m / kn) = 0.01 for the disks (m = 1), damping
    # ratio 0.5
    kn = 4e4
    run_options["osi"] = functools.partial(
        nonos.bridge.ExplicitDEMOSI, normal_stiffness=kn, normal_damping=kn**0.5
    )
run_options["Newton_max_iter"] = 1
run_options["set_external_forces"] = None
run_options["solver_options"] = options
run_options["numerics_verbose"] = False
run_options["output_contact_forces"] = False
run_options["output_frequency"] = None

with MechanicsHdf5Runner(
    mode="r+", config=runner_config, io_filename=hdf5_output
) as io:
    io.run(run_options)

# the trajectories of the two disks: rows (t, id, x, y, ...) of data/dynamic
with h5py.File(hdf5_output, "r") as f:
    dynamic = f["data/dynamic"][:]

t0 = dynamic[:, 0].min()
first = dynamic[dynamic[:, 0] == t0]
ids = [int(first[np.argmin(np.abs(first[:, 2] - p[0])), 1]) for p in positions]
a = dynamic[dynamic[:, 1] == ids[0]]
b = dynamic[dynamic[:, 1] == ids[1]]

dx = np.abs(b[:, 2] - a[:, 2] - shift).max()
dy = np.abs(b[:, 3] - a[:, 3]).max()
fall = a[0, 3] - a[:, 3].min()
print(
    f"disks {ids}: max |dx - {shift}| = {dx:.3e}, max |dy| = {dy:.3e}, fall {fall:.3e}"
)

if dx > 1e-9 or dy > 1e-9:
    print("FAILED: the two halves of the scene differ")
    sys.exit(1)
print("OK")
