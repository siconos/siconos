#!/usr/bin/env python
# siconos-tutorials GeometricPrimitives/convexhull2d_on_box.py, run with the
# vnative backend (vkernel): two convex polygons falling on a Box2d ground.
import sys

import numpy as np
import nonos.bridge as sio  # the vnative backend options
import siconos.numerics as sn
from siconos.io.mechanics_run import (
    MechanicsHdf5Runner,
    MechanicsHdf5Runner_run_options,
    RunnerConfig,
)
from siconos.mechanics.collision.convexhull import ConvexHull2d
from siconos.mechanics.collision.tools import Contactor

T = float(sys.argv[1]) if len(sys.argv) > 1 else 3.0
density = 1000.0
fname = "convexhull2d_vnative.hdf5"
config = RunnerConfig("vnative")

with MechanicsHdf5Runner(mode="w", io_filename=fname, config=config) as io:
    for name, verts, y in (
        ("ConvexHull", [[0.0, 0.0], [2.0, 0.0], [0, 2.0]], 1.0),
        ("ConvexHull2", [[0.0, 0.0], [2.0, 0.0], [0, 2.0], [1, 3.0], [3, 2.0]], 9.0),
    ):
        cm = ConvexHull2d(verts).centroid()
        verts = np.array(verts)[:] - cm[:]
        ch2d = ConvexHull2d(verts)
        inertia, area = ch2d.inertia(ch2d.centroid())
        io.add_convex_shape(name, verts)
        io.add_object(
            name.lower(),
            [Contactor(name)],
            translation=[0, y],
            orientation=[0],
            velocity=[0, 0, 0.0],
            mass=area * density,
            inertia=inertia * density,
        )
    io.add_primitive_shape(
        "Ground", "Box2d", (15, 1), insideMargin=0.0, outsideMargin=0.0
    )
    io.add_object("ground", [Contactor("Ground")], translation=[0, -0.5])
    io.add_Newton_impact_friction_nsl("contact", mu=0.1, e=0.5)

options = sn.SolverOptions(sn.solver_ids.SICONOS_FRICTION_2D_NSGS)
options.iparam[sn.params.SICONOS_IPARAM_MAX_ITER] = 1000
options.dparam[sn.params.SICONOS_DPARAM_TOL] = 1e-8

vnative_options = sio.SpaceFilterOptions()
vnative_options.neighborhood_radius = 6.0  # > largest body extent + gap
vnative_options.min_radius = 0.5  # ground segment sampling

run_options = MechanicsHdf5Runner_run_options()
run_options["t0"] = 0
run_options["T"] = T
run_options["h"] = 0.001
run_options["theta"] = 0.50001
run_options["solver_options"] = options
run_options["vnative_options"] = vnative_options
run_options["output_frequency"] = 10
run_options["output_contact_info"] = True
run_options["skip_last_update_output"] = True

with MechanicsHdf5Runner(mode="r+", io_filename=fname, config=config) as io:
    io.run(run_options)

# summary from the output
import h5py

with h5py.File(fname, "r") as f:
    dyn = f["data"]["dynamic"][:]
    last_t = dyn[:, 0].max()
    last = dyn[dyn[:, 0] == last_t]
    for row in last:
        print(
            "final t={:.3f} id={:d} x={:.4f} y={:.4f} theta={:.4f}".format(
                row[0], int(row[1]), row[2], row[3], row[4] if dyn.shape[1] > 4 else 0
            )
        )
    print("cf rows:", f["data"]["cf"].shape if "cf" in f["data"] else None)
