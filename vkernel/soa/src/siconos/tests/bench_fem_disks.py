# Benchmark: N disks stacked above a FEM square (unit square), on a ground segment.
# Models the reported slow configuration: many disks + FEM in contact.
# Timings: the io.mechanics timer lines (computeOneStep, updateInteractions...).
#
# usage: python3 bench_fem_disks.py <mesh.msh> <ndisks> [hdf5_output]
#          [neighborhood_factor] [spacing_factor] [T]
#
# neighborhood_factor: neighborhood radius of the space filter, in disk radii.
# Without it the SpaceFilterOptions defaults (radius 2.1, meant for unit
# disks) apply: with disks of radius 0.01 every pair of disks is then a
# neighbor, i.e. N(N-1)/2 interactions.
# spacing_factor: grid spacing in disk radii (default 2.05). 2.0 starts with
# every disk touching its neighbors, i.e. the whole pile in contact.
# T: final time (default 0.5, i.e. 100 steps).

import math
import sys

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

shape_filename = sys.argv[1]
ndisks = int(sys.argv[2]) if len(sys.argv) > 2 else 100
hdf5_output = sys.argv[3] if len(sys.argv) > 3 else "bench_fem_disks.hdf5"
neighborhood_factor = float(sys.argv[4]) if len(sys.argv) > 4 else None
spacing_factor = float(sys.argv[5]) if len(sys.argv) > 5 else 2.05
final_time = float(sys.argv[6]) if len(sys.argv) > 6 else 0.5

disk_radius = 0.01

backend = "vnative"
runner_config = RunnerConfig(backend)

# grid layout: rows of disks just above the unit square [0,1]^2
spacing = spacing_factor * disk_radius
nx = max(1, int(math.floor((0.9 - 2 * disk_radius) / spacing)))
ny = int(math.ceil(ndisks / nx))

x_vals = 0.05 + disk_radius + np.arange(nx) * spacing
y_vals = 1.0 + disk_radius + 1e-3 + np.arange(ny) * spacing
xx, yy = np.meshgrid(x_vals, y_vals, indexing="xy")
positions = np.column_stack([xx.ravel(), yy.ravel()])[:ndisks]
ndisks = len(positions)
print(
    f"bench: {ndisks} disks of radius {disk_radius}, grid {nx}x{ny} "
    f"(spacing {spacing:.4f}), h5 = {hdf5_output}",
    file=sys.stderr,
)

with MechanicsHdf5Runner(config=runner_config, io_filename=hdf5_output) as io:
    io.add_primitive_shape("DiskR", "Disk", [disk_radius])

    # one aggregate: the disks are created in a batch at import
    velocities = np.zeros((ndisks, 3), dtype=np.float64)
    velocities[:, 2] = -0.1
    io.add_objects(
        "disks",
        [Contactor("DiskR")],
        translations=positions,
        orientations=np.zeros(ndisks, dtype=np.float64),
        velocities=velocities,
        mass=1,
        inertia=1,
    )

    io.add_primitive_shape("Ground-1", "Segment", (-10, 0, 10, 0))

    io.add_shape_data_from_file("Square", shape_filename)
    io.add_object(
        "square",
        [Contactor("Square")],
        material=(2500, 1e11, 0.3),
    )

    io.add_Newton_impact_friction_nsl("contact", mu=0.5, e=0)

options = sn.SolverOptions(sn.solver_ids.SICONOS_FRICTION_2D_NSGS)
options.iparam[sn.params.SICONOS_IPARAM_MAX_ITER] = 100
options.dparam[sn.params.SICONOS_DPARAM_TOL] = 1e-2
# legacy FC2D error measure: the default became Euclidean scaling upstream
options.iparam[sn.params.SICONOS_FRICTION_2D_SCALING_NORM] = (
    sn.params.SICONOS_FRICTION_2D_SCALING_NO
)
options.iparam[sn.params.SICONOS_NSGS_FREEZING_CONTACT] = 10

run_options = MechanicsHdf5Runner_run_options()
run_options["with_timer"] = True
run_options["t0"] = 0
run_options["T"] = final_time
run_options["h"] = 0.005
run_options["theta"] = 0.50001
run_options["Newton_max_iter"] = 1
run_options["set_external_forces"] = None
run_options["solver_options"] = options
run_options["numerics_verbose"] = False
run_options["output_contact_forces"] = False
run_options["output_frequency"] = None

if neighborhood_factor is not None:
    vnative_options = nonos.bridge.SpaceFilterOptions()
    vnative_options.neighborhood_radius = neighborhood_factor * disk_radius
    vnative_options.min_radius = disk_radius
    run_options["vnative_options"] = vnative_options

with MechanicsHdf5Runner(
    mode="r+", config=runner_config, io_filename=hdf5_output
) as io:
    io.run(run_options)
