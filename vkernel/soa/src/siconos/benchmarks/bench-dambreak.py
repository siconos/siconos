#!/usr/bin/env python3
# Usage: bench-dambreak.py <size> [options]   (vnative backend)
#
# 2D dam break of a debris flow: a "liquid" of many small, equal,
# low-friction disks mixed with random convex rocks of several size classes
# (small and numerous to large and rare, higher friction). The material is
# dropped in a reservoir at the top of a rough slope, behind a gate. It
# settles, then the gate is removed (time of death) and the flow runs down
# the slope into a flat valley.
#
# The slope follows the USGS debris-flow flume (31 degrees, flattening at the
# bottom into a runout zone); the reservoir aspect ratio H/W is the parameter
# of the granular column collapse literature.
#
# <size> scales the whole scene: the reservoir is <size> disk diameters
# wide, the slope and the valley are proportional, so the number of bodies
# grows like <size>**2: a small size is a laptop test, a large one a
# cluster run. The counts are printed before the run.
#
# Lengths are in disk radii, as in the other benchmarks (vview draws the
# segments with a fixed thickness); --unit is the disk radius in meters and
# goes into gravity_scale, so times are physical.
#
# Example:
#   python3 bench-dambreak.py 40
#   python3 bench-dambreak.py 40 --t-settle 0.5 --T 3

import argparse
import math
import sys

import nonos.bridge as sio
import numpy as np
import siconos.numerics as sn
import siconos.simulation as sim
from siconos.mechanics.collision.tools import Contactor
from siconos.io.mechanics_run import (
    MechanicsHdf5Runner,
    MechanicsHdf5Runner_run_options,
    RunnerConfig,
    convex_hull_2d,
    polygon_mass_properties,
)

runner_config = RunnerConfig("vnative")

size = int(sys.argv[1])

p = argparse.ArgumentParser(add_help=False)
# the liquid
p.add_argument("--unit", type=float, default=0.005, help="disk radius, m")
p.add_argument("--mu-liquid", type=float, default=0.1)
# the rocks: circumdiameters as multiples of the disk diameter, and the area
# fraction of the material made of rocks. Each class gets the same total
# area, so the count of a class goes like 1 / diameter**2.
p.add_argument("--rock-sizes", type=float, nargs="+", default=[4, 8, 16, 32])
p.add_argument("--rock-fraction", type=float, default=0.4)
# distinct random shapes per class, the rocks of a class are drawn from them
p.add_argument("--rock-shapes", type=int, default=32)
p.add_argument("--mu-rock", type=float, default=0.6)
# the densities relative to water (specific gravities): the liquid disks are
# water, the rocks granite or quartz, the usual 2.65 of rock and sediment
# grains. Only their ratio matters: the masses are density * area in disk
# radii. --density sets both, as before they were separate.
p.add_argument("--density-liquid", type=float, default=1.0)
p.add_argument("--density-rock", type=float, default=2.65)
p.add_argument("--density", type=float, default=None)
# the terrain, lengths as multiples of the reservoir width
p.add_argument("--slope", type=float, default=31.0, help="degrees")
p.add_argument("--slope-length", type=float, default=8.0)
p.add_argument("--valley-length", type=float, default=6.0)
p.add_argument("--aspect-ratio", type=float, default=1.0, help="H/W")
# normal jitter of the terrain nodes, as a fraction of the segment length
p.add_argument("--roughness", type=float, default=0.25)
# time
p.add_argument("--h", type=float, default=None)
p.add_argument("--t-settle", type=float, default=None)
p.add_argument("--T", type=float, default=None)
p.add_argument("--seed", type=int, default=0)
p.add_argument("--margin", type=float, default=None)
p.add_argument("--contact-margin", type=float, default=None)
p.add_argument("--sort-frequency", type=int, default=None)
# siconos: the hdf5 tables, vtkhdf: a ParaView file instead, or both
p.add_argument(
    "--output-format", choices=("siconos", "vtkhdf", "both"), default="siconos"
)
args, _ = p.parse_known_args(sys.argv[2:])
if args.density is not None:
    args.density_liquid = args.density_rock = args.density

rng = np.random.default_rng(args.seed)
g = 9.81 / args.unit  # in disk radii per second squared

# collision groups
LIQUID, ROCK, TERRAIN = 0, 1, 2

r_d = 1.0
d_d = 2 * r_d
W = size * d_d  # reservoir width
H = args.aspect_ratio * W  # height of the settled material, about
theta = math.radians(args.slope)

# rock classes, the largest at most a quarter of the reservoir width
rock_diameters = [s * d_d for s in args.rock_sizes if s * d_d <= W / 4]
if not rock_diameters:
    sys.exit("bench-dambreak: size too small for the smallest rock class")

io_filename = "bench-dambreak-{}.hdf5".format(size)


# --- terrain --------------------------------------------------------------
# The reservoir sits on the slope between the upstream wall at x = 0 and the
# gate at x = W. The slope runs down to x = L_s, bends to horizontal over a
# transition of about W, then the valley ends on a wall.
L_s = args.slope_length * W
L_v = args.valley_length * W
y_top = L_s * math.tan(theta)
seg = rock_diameters[0]  # segment length of the rough terrain


def terrain_y(x):
    """The smooth slope, without roughness."""
    return y_top - x * math.tan(theta)


def terrain_profile():
    """Nodes of the terrain from the top of the slope to the end of the
    valley, jittered along the normal except under the reservoir. The
    terrain goes one segment past each end wall, and the walls go below it:
    disks pressed into a corner where two segments only share an endpoint
    slip out through it."""
    nodes = []
    n_slope = max(2, int(L_s / (seg * math.cos(theta))))
    for k in range(-1, n_slope + 1):
        x = L_s * k / n_slope
        nodes.append((x, terrain_y(x)))
    # transition: the angle decreases linearly to zero
    n_tr = max(4, int(W / seg))
    for k in range(1, n_tr + 1):
        a = theta * (1 - k / n_tr)
        x, y = nodes[-1]
        nodes.append((x + seg * math.cos(a), y - seg * math.sin(a)))
    x0, y0 = nodes[-1]
    n_valley = max(2, int(L_v / seg))
    for k in range(1, n_valley + 2):
        nodes.append((x0 + L_v * k / n_valley, y0))
    nodes = np.array(nodes)

    t = np.gradient(nodes, axis=0)
    t /= np.linalg.norm(t, axis=1)[:, None]
    normal = np.column_stack([-t[:, 1], t[:, 0]])
    jitter = rng.uniform(-1, 1, len(nodes)) * args.roughness * seg
    # smooth under the reservoir, so the gate and the wall close it
    jitter[nodes[:, 0] <= W + seg] = 0
    jitter[[0, -1]] = 0
    return nodes + jitter[:, None] * normal


nodes = terrain_profile()


# --- material -------------------------------------------------------------
def make_rock(diameter):
    """A random convex polygon of at most 10 vertices and circumdiameter at
    most `diameter`, centred on its area centroid: (vertices, area, polar
    second moment of area, circumradius)."""
    n = rng.integers(5, 11)
    angles = rng.uniform(0, 2 * math.pi, n)
    radii = rng.uniform(0.6, 1.0, n) * diameter / 2
    verts = convex_hull_2d(
        np.column_stack([radii * np.cos(angles), radii * np.sin(angles)])
    )
    if len(verts) < 3:
        return make_rock(diameter)
    area, centroid, polar = polygon_mass_properties(verts)
    verts = verts - centroid
    return verts, area, polar, float(np.linalg.norm(verts, axis=1).max())


material_area = W * H
area_per_class = args.rock_fraction * material_area / len(rock_diameters)

# the shapes, args.rock_shapes per class, and the rocks: indices in them
shapes = [make_rock(dia) for dia in rock_diameters for _ in range(args.rock_shapes)]
rocks = []
for c in range(len(rock_diameters)):
    a = 0.0
    while a < area_per_class:
        rocks.append(c * args.rock_shapes + rng.integers(args.rock_shapes))
        a += shapes[rocks[-1]][1]
rocks.sort(key=lambda k: -shapes[k][3])  # the largest first

# The drop zone is a rectangle of width W above the reservoir: the rocks are
# placed in it by rejection on their bounding circles, the largest first,
# then the liquid fills the free sites of a lattice, row by row. The whole
# falls into the reservoir and settles there in layers.
packing = 0.5
x_min, x_max = 0.0, W
y_min = terrain_y(0.0) + d_d
y_max = y_min + material_area / packing / W
cell = 2 * shapes[rocks[0]][3]
grid = {}


def neighbours(x, y):
    i, j = int(x // cell), int(y // cell)
    for di in (-1, 0, 1):
        for dj in (-1, 0, 1):
            yield from grid.get((i + di, j + dj), ())


def free(x, y, r):
    return all(
        (x - cx) ** 2 + (y - cy) ** 2 >= (r + cr) ** 2
        for cx, cy, cr in neighbours(x, y)
    )


def mark(x, y, r):
    grid.setdefault((int(x // cell), int(y // cell)), []).append((x, y, r))


placed_rocks = []  # (shape index, x, y)
for k in rocks:
    cr = shapes[k][3]
    for _ in range(200):
        x = rng.uniform(x_min + cr, x_max - cr)
        y = rng.uniform(y_min + cr, y_max - cr)
        if free(x, y, cr):
            break
    else:
        # no room left: on top, the drop zone grows
        x, y = rng.uniform(x_min + cr, x_max - cr), y_max + cr
        y_max = y + cr
    mark(x, y, cr)
    placed_rocks.append((k, x, y))

pitch = 1.05 * d_d
xs = np.arange(x_min + pitch / 2, x_max - r_d, pitch)
ys = np.arange(y_min + pitch / 2, y_max - r_d, pitch)
disks = np.array([(x, y) for y in ys for x in xs if free(x, y, r_d)])

y_wall_top = y_max + 2 * d_d

print(
    "bench-dambreak: size={} disks={} rocks={} classes={} reservoir={:.0f}r "
    "slope={:.0f}r valley={:.0f}r terrain segments={} "
    "densities liquid={:g} rock={:g}".format(
        size,
        len(disks),
        len(placed_rocks),
        len(rock_diameters),
        W,
        L_s,
        L_v,
        len(nodes) - 1,
        args.density_liquid,
        args.density_rock,
    ),
    flush=True,
)

# --- time -----------------------------------------------------------------
# A body must not move more than half a disk radius in a step, or the disks
# tunnel through the terrain segments: the bound is the frictionless speed
# at the valley floor, reached from the top of the drop zone.
v_max = math.sqrt(2 * g * (y_max - nodes[:, 1].min()))
h = args.h if args.h is not None else 0.5 * r_d / v_max
# settle: twice the free fall from the top of the drop zone; flow: twice the
# frictionless run down the slope from rest
t_fall = math.sqrt(2 * (y_max - terrain_y(W)) / g)
t_settle = args.t_settle if args.t_settle is not None else 2 * t_fall
t_flow = 2 * math.sqrt(2 * L_s / math.cos(theta) / (g * math.sin(theta)))
T = args.T if args.T is not None else t_settle + t_flow

# --- scene ----------------------------------------------------------------
with MechanicsHdf5Runner(io_filename=io_filename, config=runner_config) as io:
    statics = {
        "terrain-{}".format(k): (*a, *b)
        for k, (a, b) in enumerate(zip(nodes, nodes[1:]))
    }
    # walls and gate start below the terrain, see terrain_profile
    x_end, y_valley = nodes[-2][0], nodes[-1][1]
    statics["wall-up"] = (0, terrain_y(0) - seg, 0, y_wall_top)
    statics["wall-end"] = (x_end, y_valley - seg, x_end, y_valley + H)
    statics["gate"] = (W, terrain_y(W) - seg, W, y_wall_top)
    for name, s in statics.items():
        io.add_primitive_shape(name.capitalize(), "Segment", s)
        io.add_object(
            name,
            [Contactor(name.capitalize(), collision_group=TERRAIN)],
            translation=[0, 0],
            # the dam breaks
            time_of_death=t_settle if name == "gate" else -1,
        )

    io.add_primitive_shape("Liquid", "Disk", [r_d], insideMargin=0.0, outsideMargin=0.0)
    io.add_objects(
        "liquid",
        Contactor("Liquid", collision_group=LIQUID),
        disks,
        density=args.density_liquid,
    )

    for k, (verts, _, _, _) in enumerate(shapes):
        io.add_convex_shape("Rock-{}".format(k), verts, insideMargin=0.001 * r_d)
    placed = np.array(placed_rocks)
    io.add_objects(
        "rocks",
        [
            Contactor("Rock-{}".format(k), collision_group=ROCK)
            for k in range(len(shapes))
        ],
        placed[:, 1:],
        orientations=rng.uniform(0, 2 * math.pi, len(placed)),
        shape_indices=placed[:, 0].astype(int),
        density=args.density_rock,
    )

    # the liquid lubricates every contact it is in
    for name, mu, g1, g2 in (
        ("liquid-liquid", args.mu_liquid, LIQUID, LIQUID),
        ("liquid-rock", args.mu_liquid, LIQUID, ROCK),
        ("liquid-terrain", args.mu_liquid, LIQUID, TERRAIN),
        ("rock-rock", args.mu_rock, ROCK, ROCK),
        ("rock-terrain", args.mu_rock, ROCK, TERRAIN),
    ):
        io.add_Newton_impact_friction_nsl(
            name, mu=mu, e=0, collision_group1=g1, collision_group2=g2
        )

# --- run ------------------------------------------------------------------
options = sn.SolverOptions(sn.solver_ids.SICONOS_FRICTION_2D_NSGS)
options.iparam[sn.params.SICONOS_IPARAM_MAX_ITER] = 100
options.dparam[sn.params.SICONOS_DPARAM_TOL] = 1e-2
options.iparam[sn.params.SICONOS_NSGS_FREEZING_CONTACT] = 10

run_options = MechanicsHdf5Runner_run_options()
run_options["t0"] = 0
run_options["T"] = T
run_options["h"] = h
run_options["gravity_scale"] = args.unit
run_options["solver_options"] = options
run_options["Newton_options"] = sim.LINEAR
run_options["skip_last_update_output"] = True
run_options["skip_reset_lambdas"] = True
run_options["violation_verbose"] = True
# about 100 frames per second of simulated time
run_options["output_frequency"] = max(1, int(round(0.01 / h)))
run_options["with_timer"] = True
run_options["time_stepping"] = None
run_options["output_contact_info"] = True
run_options["output_contact_work"] = False
run_options["output_format"] = args.output_format

vnative_options = sio.SpaceFilterOptions()
# The broadphase searches body centres with one radius: two bodies can touch
# when their centres are closer than the sum of their circumradii, so the
# largest rock sets it for every body.
margin = args.margin if args.margin is not None else r_d
vnative_options.neighborhood_radius = 2 * shapes[rocks[0]][3] + margin
# A disk then finds hundreds of others: only the pairs whose bounding circles
# are closer than the contact margin get an interaction (<= 0: all of them).
contact_margin = args.contact_margin if args.contact_margin is not None else margin
if contact_margin > 0:
    vnative_options.contact_margin = contact_margin
vnative_options.min_radius = r_d
if args.sort_frequency is not None:
    vnative_options.sort_frequency = args.sort_frequency
vnative_options.run_vnative_options = True
run_options["vnative_options"] = vnative_options

print(
    "bench-dambreak: h={:.2e} t_settle={:.3f} T={:.3f} steps={}".format(
        h, t_settle, T, int(round(T / h))
    ),
    flush=True,
)

with MechanicsHdf5Runner(
    io_filename=io_filename, mode="r+", config=runner_config
) as io:
    io.run(run_options)
