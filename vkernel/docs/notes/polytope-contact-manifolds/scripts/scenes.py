#!/usr/bin/env python
# Polytope contact scenarios run with the vnative backend (Moreau-Jean),
# every step written to <scene>.hdf5 with the active contacts (cf_data).
#
#   python3 scenes.py <scene> [T]
#
# scenes: flat, tilted, overhang, tipping, junction, stack, corner
import math
import sys

import numpy as np
import nonos.bridge as sio
import siconos.numerics as sn
from siconos.io.mechanics_run import (
    MechanicsHdf5Runner,
    MechanicsHdf5Runner_run_options,
    RunnerConfig,
)
from siconos.mechanics.collision.convexhull import ConvexHull2d
from siconos.mechanics.collision.tools import Contactor

scene = sys.argv[1]
density = 1000.0


def rect(w, hgt):
    return [[-w / 2, -hgt / 2], [w / 2, -hgt / 2], [w / 2, hgt / 2], [-w / 2, hgt / 2]]


def add_poly(io, name, verts, pos, theta=0.0, vel=(0.0, 0.0, 0.0)):
    verts = np.array(verts, dtype=float)
    cm = ConvexHull2d(verts).centroid()
    verts = verts - cm
    ch = ConvexHull2d(verts)
    inertia, area = ch.inertia(ch.centroid())
    io.add_convex_shape(name + "_shape", verts)
    io.add_object(
        name,
        [Contactor(name + "_shape")],
        translation=list(pos),
        orientation=[theta],
        velocity=list(vel),
        mass=area * density,
        inertia=inertia * density,
    )


def add_segment(io, name, x1, y1, x2, y2):
    io.add_primitive_shape(name + "_shape", "Segment", (x1, y1, x2, y2))
    io.add_object(name, [Contactor(name + "_shape")], translation=[0, 0])


def lowest_corner_y(w, hgt, theta):
    c, s = math.cos(theta), math.sin(theta)
    return min(s * x + c * y for x, y in rect(w, hgt))


# box of the scenes: 1.0 x 0.4
W, H = 1.0, 0.4
T = 1.0
e = 0.0
mu = 0.3

fname = f"{scene}.hdf5"
config = RunnerConfig("vnative")
with MechanicsHdf5Runner(mode="w", io_filename=fname, config=config) as io:
    if scene == "flat":
        # flat landing: both manifold points reach the ground together
        e = 0.5
        add_segment(io, "ground", -3, 0, 3, 0)
        add_poly(io, "box", rect(W, H), (0.0, H / 2 + 0.05))
    elif scene == "tilted":
        # corner-first landing: one point active first, then the second
        e = 0.5
        T = 1.5
        th = math.radians(20)
        add_segment(io, "ground", -3, 0, 3, 0)
        add_poly(io, "box", rect(W, H), (0.0, -lowest_corner_y(W, H, th) + 0.05), th)
    elif scene == "overhang":
        # resting over the end of a segment, centre of mass above it
        add_segment(io, "ground", -3, 0, 0, 0)
        add_poly(io, "box", rect(W, H), (-0.3, H / 2 + 0.01))
    elif scene == "tipping":
        # centre of mass beyond the end of the segment: the box tips over it
        T = 0.8
        add_segment(io, "ground", -3, 0, 0, 0)
        add_poly(io, "box", rect(W, H), (0.15, H / 2 + 0.01))
    elif scene == "junction":
        # two collinear segments (like two mesh edges): one pair per segment
        add_segment(io, "ground1", -3, 0, 0, 0)
        add_segment(io, "ground2", 0, 0, 3, 0)
        add_poly(io, "box", rect(W, H), (0.1, H / 2 + 0.01))
    elif scene == "stack":
        # polytope on polytope: a small box on the edge of a large one
        add_segment(io, "ground", -3, 0, 3, 0)
        add_poly(io, "base", rect(2.0, 0.4), (0.0, 0.2))
        add_poly(io, "top", rect(0.6, 0.3), (0.85, 0.4 + 0.15 + 0.02))
    elif scene == "corner":
        # a square falling on its corner onto a large box (vertex on face)
        T = 1.5
        add_segment(io, "ground", -3, 0, 3, 0)
        add_poly(io, "base", rect(2.0, 0.4), (0.0, 0.2))
        th = math.radians(30)
        add_poly(
            io,
            "square",
            rect(0.4, 0.4),
            (0.0, 0.4 - lowest_corner_y(0.4, 0.4, th) + 0.05),
            th,
        )
    else:
        raise SystemExit(f"unknown scene {scene}")
    io.add_Newton_impact_friction_nsl("contact", mu=mu, e=e)

T = float(sys.argv[2]) if len(sys.argv) > 2 else T

options = sn.SolverOptions(sn.solver_ids.SICONOS_FRICTION_2D_NSGS)
options.iparam[sn.params.SICONOS_IPARAM_MAX_ITER] = 10000
options.dparam[sn.params.SICONOS_DPARAM_TOL] = 1e-12

vnative_options = sio.SpaceFilterOptions()
vnative_options.neighborhood_radius = 3.0
vnative_options.min_radius = 0.1  # segment sampling

run_options = MechanicsHdf5Runner_run_options()
run_options["t0"] = 0
run_options["T"] = T
run_options["h"] = 1e-3
run_options["theta"] = 0.5
run_options["constraint_activation_threshold"] = 0.0
run_options["solver_options"] = options
run_options["vnative_options"] = vnative_options
run_options["output_frequency"] = 1
run_options["output_contact_info"] = True
run_options["skip_last_update_output"] = True
run_options["verbose"] = False

with MechanicsHdf5Runner(mode="r+", io_filename=fname, config=config) as io:
    io.run(run_options)
