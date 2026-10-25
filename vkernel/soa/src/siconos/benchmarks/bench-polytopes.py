#!/usr/bin/env python3
# Usage: bench-polytopes.py <backend> <N> [options]
#
# N x N grid of random small convex-polygon grains in a box under gravity,
# run through a backend (vnative or bullet).
#
# Modes:
#   default            : settle the grains (full sim), original behaviour.
#   --dump PATH        : after the run, dump the final body positions + shape
#                        vertices to PATH (npz). Use the vnative run to produce a
#                        settled static configuration.
#   --load PATH --compare-steps K : rebuild the scene from a dumped PATH,
#                        place bodies at the dumped positions with zero velocity,
#                        and run K time-steps only. Reports per-step timing so
#                        vnative and bullet can be compared on the SAME geometry.
#
# Example (settle, then compare):
#   python3 bench-polytopes.py vnative 30 --dump settled-30.npz
#   python3 bench-polytopes.py bullet  30 --load settled-30.npz --compare-steps 20
#   python3 bench-polytopes.py vnative  30 --load settled-30.npz --compare-steps 20

import argparse
import math
import random
import sys
import time

import numpy as np
import siconos.numerics as sn
import siconos.simulation as sim
from siconos.mechanics.collision.convexhull import ConvexHull2d
from siconos.mechanics.collision.tools import Contactor
from siconos.io.mechanics_run import (
    MechanicsHdf5Runner,
    MechanicsHdf5Runner_run_options,
    RunnerConfig,
)

# ---------------------------------------------------------------------------

backend = str(sys.argv[1])
runner_config = RunnerConfig(backend)

N = int(sys.argv[2])

# Parse our own optional flags (ignored by siconos internals).
rest = sys.argv[3:]
p = argparse.ArgumentParser(add_help=False)
p.add_argument("--dump", metavar="PATH")
p.add_argument("--load", metavar="PATH")
p.add_argument("--compare-steps", type=int, default=0)
p.add_argument("--h", type=float, default=0.005)
# neighborhood radius = 2 * largest grain circumradius + margin; the margin
# plays the role of bullet's contactBreakingThreshold below
p.add_argument("--margin", type=float, default=0.4)
p.add_argument("--neighborhood-radius", type=float, default=None)
args, _ = p.parse_known_args(rest)

random.seed(0)
np.random.seed(0)

density = 1000.0
grain_size = 1.0

io_filename = "bench-polytopes-{}-{}.hdf5".format(backend, N)
# Compare runs (load + timed steps) write to a separate file so they never
# append to / clobber the long settle run's hdf5.
if args.load or args.compare_steps:
    io_filename = "bench-polytopes-{}-{}-cmp.hdf5".format(backend, N)

base = 2 * N
dbase = base
offset = base * 0.0
vinit = base * 0.0
h1 = 4
h2 = 7
l1 = 9
l2 = 11


def make_grain_verts():
    """A random (but repeatable) 4-vertex convex polygon, centred on its
    centroid, scaled to grain_size."""
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
    verts = np.array(vertices)
    scale = grain_size / max(np.ptp(verts, axis=0))
    cm = ConvexHull2d(verts).centroid()
    verts = (verts - cm) * scale
    return verts, ConvexHull2d(verts)


# Loaded state (set when --load is given). verts_source[k] -> (verts, ch).
q_dump = None
verts_source = None
circumradii = []


def build_scene(io, place_func):
    """Build the box + N grains. place_func(idx, i, j, cname, verts, ch) adds
    each grain, using either the freshly generated or the dumped geometry."""
    io.add_primitive_shape(
        "Ground-1", "Segment", (-base, h1 * base + dbase, 0, h1 * base)
    )
    io.add_primitive_shape("Ground-2", "Segment", (0, h1 * base, l1 * base, 0))
    io.add_primitive_shape("Ground-3", "Segment", (l1 * base, 0, l2 * base, 0))
    io.add_primitive_shape("Wall-1", "Segment", (l2 * base, 0, l2 * base, h2 * base))
    io.add_primitive_shape(
        "Wall-2", "Segment", (-base, h1 * base + dbase, -base, h2 * base)
    )

    io.add_object("ground-1", [Contactor("Ground-1")], translation=[0, 0])
    io.add_object("ground-2", [Contactor("Ground-2")], translation=[0, 0])
    io.add_object("ground-3", [Contactor("Ground-3")], translation=[0, 0])
    io.add_object("wall-1", [Contactor("Wall-1")], translation=[0, 0])
    io.add_object("wall-2", [Contactor("Wall-2")], translation=[0, 0])

    count = 0
    n_side = int(math.ceil(math.sqrt(N)))
    for i in range(n_side):
        for j in range(n_side):
            if count >= N:
                return count
            if verts_source is None:
                verts, ch = make_grain_verts()
            else:
                verts, ch = verts_source[count]
            circumradii.append(float(np.linalg.norm(verts, axis=1).max()))
            cname = "PolyCS_{}_{}".format(i, j)
            io.add_convex_shape(cname, verts, insideMargin=0.001 * grain_size)
            place_func(count, i, j, cname, verts, ch)
            count += 1
        if count >= N:
            return count
    return count


def make_place_grid():
    def place(idx, i, j, cname, verts, ch):
        inertia, area = ch.inertia(ch.centroid())
        io.add_object(
            "poly{}-{}".format(i, j),
            [Contactor(cname)],
            translation=[-base + 2 * i + 1, h1 * base + 2 * j + 1 + dbase + offset],
            orientation=[0],
            velocity=[vinit, 0, 0],
            mass=area * density,
            inertia=inertia * density,
        )

    return place


def make_place_loaded():
    def place(idx, i, j, cname, verts, ch):
        dq = q_dump[idx]
        inertia, area = ch.inertia(ch.centroid())
        io.add_object(
            "poly{}-{}".format(i, j),
            [Contactor(cname)],
            translation=[float(dq[0]), float(dq[1])],
            orientation=[float(dq[2])],
            velocity=[0, 0, 0],
            mass=area * density,
            inertia=inertia * density,
        )

    return place


# --- load phase: read dumped positions + vertices -----------------------
if args.load:
    blob = np.load(args.load, allow_pickle=True)
    q_dump = blob["q"]
    raw_verts = blob["verts"]
    verts_source = [
        (
            np.asarray(raw_verts[k], dtype=float),
            ConvexHull2d(np.asarray(raw_verts[k], dtype=float)),
        )
        for k in range(len(q_dump))
    ]
    if len(q_dump) != N:
        sys.exit(
            "bench-polytopes: --load has {} bodies but N={} requested".format(
                len(q_dump), N
            )
        )


with MechanicsHdf5Runner(io_filename=io_filename, config=runner_config) as io:
    if args.load:
        build_scene(io, make_place_loaded())
    else:
        build_scene(io, make_place_grid())

    io.add_Newton_impact_friction_nsl("contact", mu=0.5, e=0)


options = sn.SolverOptions(sn.solver_ids.SICONOS_FRICTION_2D_NSGS)
options.iparam[sn.params.SICONOS_IPARAM_MAX_ITER] = 100
options.dparam[sn.params.SICONOS_DPARAM_TOL] = 1e-2
options.iparam[sn.params.SICONOS_NSGS_FREEZING_CONTACT] = 10

run_options = MechanicsHdf5Runner_run_options()
run_options["t0"] = 0
# without it the runner keeps its default h = 0.0005 and T = K * args.h runs
# 10 K steps
run_options["h"] = args.h
# Stabilization gravity. The original density-scaling gravity_scale=1/N is far
# too weak for the tall ramp scene: at N=30 gravity ~0.33 and grains starting
# at y~300 never reach the floor within T=10, so the dumped "rest" config is
# still in free fall and the comparison times near-vacuum motion instead of
# contact resolution. Use full gravity for the settle run so grains come to
# rest on the ramp, then dump that genuinely static configuration.
if args.compare_steps > 0:
    run_options["T"] = args.compare_steps * args.h
    run_options["gravity_scale"] = 1.0
else:
    run_options["T"] = 10
    run_options["gravity_scale"] = 1.0

run_options["solver_options"] = options
run_options["constraint_activation_threshold"] = 1e-05

# Bullet backend: configure collision options for convex polytopes (without
# this the default gives missed contacts and blocks fall through the ramp).
# NOTE: bullet_options MUST be instantiated before its fields are set; an earlier
# version referenced bullet_options before assignment, raising NameError, which
# the bare except swallowed so the run fell back to SiconosBulletOptions defaults.
try:
    import siconos.mechanics.collision.bullet as sb

    bullet_options = sb.SiconosBulletOptions()
    bullet_options.worldScale = 0.01
    bullet_options.collisionMargin = 1.0
    bullet_options.contactBreakingThreshold = 0.4
    bullet_options.dimension = sb.TwoD
    bullet_options.perturbationIterations = 3
    bullet_options.minimumPointsPerturbationThreshold = 3
    run_options["bullet_options"] = bullet_options
except Exception:
    pass
run_options["Newton_options"] = sim.LINEAR

run_options["skip_last_update_output"] = True
run_options["skip_reset_lambdas"] = True

run_options["violation_verbose"] = True
run_options["output_frequency"] = 1
run_options["with_timer"] = True
run_options["time_stepping"] = None
run_options["output_contact_info"] = True
run_options["output_contact_work"] = False

import nonos.bridge as sio  # noqa: E402

# per-phase timing of the vnative space filter (search, update_index_set0,
# narrowphase, sort): seconds per phase, one entry per step
ui_phases = {}


def timed_update_interactions(self, step):
    if not self._initialized:
        ui_phases.clear()
        return untimed_update_interactions(self, step)
    t = [time.perf_counter()]
    self._ngbh.update(step)
    t.append(time.perf_counter())
    self._ngbh.search()
    t.append(time.perf_counter())
    self._handle.update_index_set0(step)
    t.append(time.perf_counter())
    self._narrowphase.ensure(step)
    t.append(time.perf_counter())
    freq = self._options.sort_frequency
    if freq > 0 and step % freq == 0:
        self._ngbh.sort()
    t.append(time.perf_counter())
    for name, a, b in zip(
        ("update", "search", "update_index_set0", "narrowphase", "sort"), t, t[1:]
    ):
        ui_phases.setdefault(name, []).append(b - a)


untimed_update_interactions = sio.SpaceFilter.updateInteractions
sio.SpaceFilter.updateInteractions = timed_update_interactions

vnative_options = sio.SpaceFilterOptions()
p_sort = argparse.ArgumentParser(add_help=False)
p_sort.add_argument("--sort-frequency", type=int, default=None)
args_sort, _ = p_sort.parse_known_args(rest)
if args_sort.sort_frequency is not None:
    vnative_options.sort_frequency = args_sort.sort_frequency
# The broadphase searches grain centres only: two grains can touch when
# their centres are closer than the sum of their circumradii. A larger
# radius only adds candidate pairs, each bound to two interactions that the
# narrowphase and every compute_h / compute_jachq then process.
if args.neighborhood_radius is not None:
    vnative_options.neighborhood_radius = args.neighborhood_radius
else:
    vnative_options.neighborhood_radius = 2 * max(circumradii) + args.margin
print(
    "vnative neighborhood_radius = {:.3f}".format(vnative_options.neighborhood_radius)
)
vnative_options.min_radius = 1.0
vnative_options.run_vnative_options = True

run_options["vnative_options"] = vnative_options

if args.compare_steps > 0 and args.load is None:
    sys.exit("bench-polytopes: --compare-steps requires --load PATH")

if args.compare_steps > 0:
    print(
        "=== compare: backend={} N={} steps={} ===".format(
            backend, N, args.compare_steps
        ),
        flush=True,
    )
    # the hdf5 output is not part of the step: write only the first steps
    run_options["output_frequency"] = 0
    t0 = time.perf_counter()
    with MechanicsHdf5Runner(
        io_filename=io_filename, mode="r+", config=runner_config
    ) as io:
        io.run(run_options)
        timing = dict(io._timing)
    elapsed = time.perf_counter() - t0

    # per-step means over the timed loop only (setup excluded). vnative runs
    # the space filter + narrowphase in updateInteractions, bullet inside
    # computeOneStep: their sum is the comparable step cost.
    def mean_ms(name):
        v = timing.get(name, [])
        return 1000.0 * sum(v) / len(v) if v else 0.0

    nsteps = len(timing.get("computeOneStep", []))
    ui_ms = mean_ms("updateInteractions")
    cos_ms = mean_ms("computeOneStep")
    print(
        "compare: elapsed={:.4f}s steps={} updateInteractions={:.3f}ms "
        "computeOneStep={:.3f}ms per_step={:.3f}ms".format(
            elapsed, nsteps, ui_ms, cos_ms, ui_ms + cos_ms
        ),
        flush=True,
    )
    if ui_phases:
        print(
            "compare: space filter "
            + " ".join(
                "{}={:.3f}ms".format(k, 1000.0 * sum(v) / len(v))
                for k, v in ui_phases.items()
            ),
            flush=True,
        )
    sys.exit(0)

with MechanicsHdf5Runner(
    io_filename=io_filename, mode="r+", config=runner_config
) as io:
    io.run(run_options)

# Dump final positions + shape vertices for cross-backend comparison.
if args.dump:
    with MechanicsHdf5Runner(
        io_filename=io_filename, mode="r", config=runner_config
    ) as io:
        dyn = io.dynamic_data()[:]
    # rows: [time, ds_id, x, y, z, qw, qx, qy, qz]; keep the last complete
    # step for each body, ordered by ds id (1..N). 2D q is (x, y, theta),
    # theta = 2 atan2(qz, qw) (column 4 is z, always 0 in 2D).
    if dyn.shape[0] > 0:
        last_time = dyn[:, 0].max()
        mask = dyn[:, 0] == last_time
        rows = dyn[mask]
        order = np.argsort(rows[:, 1])
        rows = rows[order]
        ids = rows[:, 1].astype(int)
        theta = 2.0 * np.arctan2(rows[:, 8], rows[:, 5])
        qs = np.column_stack([rows[:, 2], rows[:, 3], theta])
    else:
        ids = np.array([], dtype=int)
        qs = np.empty((0, 3))

    # Re-collect the per-grain shape vertices in construction order. The RNG
    # sequence is deterministic, so replaying make_grain_verts() reproduces it.
    dumped_verts = np.empty(N, dtype=object)
    random.seed(0)
    np.random.seed(0)
    for k in range(N):
        v, _ = make_grain_verts()
        dumped_verts[k] = v

    np.savez(
        args.dump,
        q=qs,
        verts=dumped_verts,
        ids=ids,
        n=N,
        backend=backend,
    )
    print(
        "=== dump: backend={} N={} bodies={} -> {} ===".format(
            backend, N, len(ids), args.dump
        ),
        flush=True,
    )
