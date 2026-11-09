#!/usr/bin/env python
# siconos-tutorials GranularMaterial/2d_flow.py with the vnative backend
# (vkernel): random convex quadrilateral grains falling on a Box2d ground.
# The grains are all created at t = 0 (no time_of_birth).
#
# usage: 2d_flow_vnative.py [n_row n_col T]
import math
import random
import sys

import nonos.bridge as sio  # the vnative backend options
import numpy
import siconos.numerics as sn
from siconos.io.mechanics_run import (
    MechanicsHdf5Runner,
    MechanicsHdf5Runner_run_options,
    RunnerConfig,
)
from siconos.mechanics.collision.convexhull import ConvexHull2d
from siconos.mechanics.collision.tools import Contactor

n_row = int(sys.argv[1]) if len(sys.argv) > 1 else 5
n_col = int(sys.argv[2]) if len(sys.argv) > 2 else 3
T = float(sys.argv[3]) if len(sys.argv) > 3 else 1.0
hstep = 1e-3
density = 1000.0
grain_size = 0.1
fname = "2d_flow_vnative.hdf5"
config = RunnerConfig("vnative")
random.seed(0)
numpy.random.seed(0)


def create_grain(io, name, cname, grain_size, density, trans):
    rd = [math.pi / 2 * random.gauss(0.5, 0.2) for _ in range(16)]

    def vert(id1, id2, a, b):
        return (
            a * math.cos(rd[id1]) * math.cos(rd[id2]),
            b * math.sin(rd[id1]) * math.cos(rd[id2]),
        )

    vertices = [
        vert(0, 1, 1, 1),
        vert(2, 3, 1, -1),
        vert(4, 5, -1, 1),
        vert(6, 7, -1, -1),
    ]
    scale = grain_size / max(
        numpy.array(vertices).max(axis=0) - numpy.array(vertices).min(axis=0)
    )
    cm = ConvexHull2d(vertices).centroid()
    vertices = (numpy.array(vertices)[:] - cm[:]) * scale
    ch2d = ConvexHull2d(vertices)
    io.add_convex_shape(cname, vertices, insideMargin=0.001 * grain_size)
    inertia, area = ch2d.inertia(ch2d.centroid())
    io.add_object(
        name,
        [Contactor(cname)],
        translation=trans,
        orientation=[0],
        velocity=[0, 0, 0],
        mass=area * density,
        inertia=inertia * density,
    )


with MechanicsHdf5Runner(mode="w", io_filename=fname, config=config) as io:
    sizes = numpy.random.uniform(
        low=grain_size - 0.1 / 2, high=grain_size + 0.1 / 2, size=n_row * n_col
    )
    k = 0
    x_shift = 2.0
    for i in range(n_row):
        for j in range(n_col):
            trans = [(i - n_row / 2.0) * x_shift * grain_size, j * x_shift * grain_size]
            create_grain(
                io,
                "rock_{}_{}".format(i, j),
                "RockCS_{}_{}".format(i, j),
                sizes[k],
                density,
                trans,
            )
            k += 1
    io.add_primitive_shape(
        "Ground", "Box2d", (25, 0.1), insideMargin=0.0, outsideMargin=0.0
    )
    io.add_object("ground", [Contactor("Ground")], translation=[0, -0.5 - 0.2])
    io.add_Newton_impact_friction_nsl("contact", mu=0.1, e=0.5)

options = sn.SolverOptions(sn.solver_ids.SICONOS_FRICTION_2D_NSGS)
options.iparam[sn.params.SICONOS_IPARAM_MAX_ITER] = 1000
options.dparam[sn.params.SICONOS_DPARAM_TOL] = 1e-8

vnative_options = sio.SpaceFilterOptions()
vnative_options.neighborhood_radius = 0.4  # > 2 x largest grain extent (0.15)
vnative_options.min_radius = 0.05  # ground segment sampling

run_options = MechanicsHdf5Runner_run_options()
run_options["t0"] = 0
run_options["T"] = T
run_options["h"] = hstep
run_options["theta"] = 0.50001
run_options["solver_options"] = options
run_options["vnative_options"] = vnative_options
run_options["output_frequency"] = 100
run_options["output_contact_info"] = True
run_options["skip_last_update_output"] = True
run_options["with_timer"] = True

with MechanicsHdf5Runner(mode="r+", io_filename=fname, config=config) as io:
    io.run(run_options)

import h5py  # noqa: E402

with h5py.File(fname, "r") as f:
    dyn = f["data"]["dynamic"][:]
    last = dyn[dyn[:, 0] == dyn[:, 0].max()]
    print(
        "final t={:.3f} grains={} y min/mean/max = {:.4f} {:.4f} {:.4f}, "
        "|v| check skipped".format(
            last[0, 0], len(last), last[:, 3].min(), last[:, 3].mean(), last[:, 3].max()
        )
    )
