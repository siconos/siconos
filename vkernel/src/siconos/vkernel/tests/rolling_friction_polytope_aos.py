# Siconos is a program dedicated to modeling, simulation and control
# of non smooth dynamical systems.
#
# SPDX-FileCopyrightText: INRIA
# SPDX-License-Identifier: Apache-2.0
#

# Check: rolling resistance on polytopes, vnative backend against the
# AOS kernel. The reference is polygon_line_rolling_aos: a
# regular polygon (m = 1, g = 9.81) on a line, one Lagrangian2d3DR
# interaction per vertex (normal, tangent, rolling rows), the slope given by
# tilting the gravity. As in the AOS kernel, each contact point of a
# polytope resists rolling: |m_r| <= mu_r r_n. The vnative run puts the
# polygon on a tilted segment and its trajectory is read in the frame of the
# segment. The cases:
#
#   slide        a square sliding on a face: it stays flat and stops at
#                v0^2 / (2 mu g), the rolling rows do not change it
#   tip-<mu_r>   a square launched pivoting on a corner, above the tipping
#                threshold without rolling resistance: it tips over with
#                mu_r = 0, falls back with mu_r = 0.1
#   octagon-<mu_r> an octagon rolling over its vertices
#   incline      a square at rest on a slope
#
# The hdf5 tables are float32: the tolerance is on that scale. Each vnative
# run is in its own process (a second simulation in the same process
# crashes the vnative backend).
#
# usage: python3 rolling_friction_polytope_aos.py <polygon_line_rolling_aos>
# exit status 0 if every case agrees

import math
import os
import sys
import tempfile

import h5py
import numpy as np

# the two kernels are compared at the same h, coarse where they agree at
# any h. On a face (slide, incline) the two kernels may stop the square one
# step apart, an error of order g h^2, and a coarse h rocks it: h = 1e-3
# there. An output every OUTPUT_DT (s): the writing of an output costs more
# than a step.
H = 5e-3
H_FACE = 1e-3
OUTPUT_DT = 0.02
m, g = 1.0, 9.81
SQUARE = 1 / math.sqrt(2)  # the circumradius of the unit square
# above the tipping threshold of the unit square pivoting on a corner,
# omega^2 >= 2 m g (Rc - 1/2) / (I + m Rc^2) = 6.09
W_TIP = 2.6

# the cases: (N, Rc, mu, mu_r, tan alpha, vx0, vy0, w0, T), the initial
# velocity in the frame of the line
CASES = {
    "slide": (4, SQUARE, 0.3, 0.1, 0.0, 2.0, 0.0, 0.0, 0.8),
    "tip-0": (4, SQUARE, 0.5, 0.0, 0.0, W_TIP / 2, W_TIP / 2, -W_TIP, 1.0),
    "tip-0.1": (4, SQUARE, 0.5, 0.1, 0.0, W_TIP / 2, W_TIP / 2, -W_TIP, 1.0),
    "octagon-0": (8, 1.0, 0.5, 0.0, 0.0, 3.0, 0.0, -3.0 / math.cos(math.pi / 8), 2.0),
    "octagon-0.05": (
        8,
        1.0,
        0.5,
        0.05,
        0.0,
        3.0,
        0.0,
        -3.0 / math.cos(math.pi / 8),
        2.0,
    ),
    "incline": (4, SQUARE, 0.5, 0.1, 0.2, 0.0, 0.0, 0.0, 0.5),
}


def step(case):
    """the h of a case"""
    return H_FACE if case in ("slide", "incline") else H


def vertices(N, Rc):
    """the vertices in the frame of the polygon, as the reference: the
    first and the last on the line at t = 0, mirror images of each other"""
    v = np.zeros((N, 2))
    for k in range((N + 1) // 2):
        phi = -math.pi / 2 + math.pi / N + 2 * math.pi * k / N
        v[k] = Rc * math.cos(phi), Rc * math.sin(phi)
        v[N - 1 - k] = -v[k, 0], v[k, 1]
    return v


def run_vnative(case, filename):
    import siconos.numerics as sn
    from siconos.io.mechanics_run import (
        MechanicsHdf5Runner,
        MechanicsHdf5Runner_run_options,
        RunnerConfig,
    )
    from siconos.mechanics.collision.tools import Contactor

    import nonos  # noqa: F401 registers the vnative backend

    N, Rc, mu, mu_r, tan, vx0, vy0, w0, T = CASES[case]
    alpha = math.atan(tan)
    c, s = math.cos(alpha), math.sin(alpha)
    inertia = m * Rc**2 / 6 * (1 + 2 * math.cos(math.pi / N) ** 2)
    # the reference frame tilted by -alpha: the slope goes down along +x
    L = 50.0
    d = -vertices(N, Rc)[0, 1]
    # exactly on a tilted segment, the round-off of the gap (the input tables
    # are float32) can leave a contact inactive for one step: the polygon
    # starts 1e-7 inside it
    if tan:
        d -= 1e-7
    config = RunnerConfig("vnative")
    with MechanicsHdf5Runner(config=config, io_filename=filename) as io:
        io.add_convex_shape("Polygon", vertices(N, Rc))
        io.add_primitive_shape("Ground", "Segment", (-L * c, L * s, L * c, -L * s))
        io.add_Newton_impact_rolling_friction_nsl("contact", mu=mu, mu_r=mu_r, e=0.0)
        io.add_object(
            "polygon",
            [Contactor("Polygon")],
            translation=[d * s, d * c],
            orientation=[-alpha],
            velocity=[vx0 * c + vy0 * s, -vx0 * s + vy0 * c, w0],
            mass=m,
            inertia=inertia,
        )
        io.add_object("ground", [Contactor("Ground")], translation=[0, 0])

    options = sn.SolverOptions(sn.solver_ids.SICONOS_ROLLING_FRICTION_2D_NSGS)
    options.iparam[sn.params.SICONOS_IPARAM_MAX_ITER] = 1000
    options.dparam[sn.params.SICONOS_DPARAM_TOL] = 1e-12
    run_options = MechanicsHdf5Runner_run_options()
    run_options["t0"] = 0
    run_options["T"] = T
    run_options["h"] = step(case)
    run_options["theta"] = 0.5
    run_options["solver_options"] = options
    run_options["output_frequency"] = round(OUTPUT_DT / step(case))
    # the runner divides its gravity, 9.80665, by gravity_scale
    run_options["gravity_scale"] = 9.80665 / g
    with MechanicsHdf5Runner(mode="r+", config=config, io_filename=filename) as io:
        io.run(run_options)


def vnative(case, filename):
    """t, then x, y, theta, vx, vy, w in the frame of the line, from the file
    of the vnative run"""
    with h5py.File(filename, "r") as f:
        q, v = f["data/dynamic"][:], f["data/velocities"][:]
    ds = v[0, 1]
    q, v = q[q[:, 1] == ds], v[v[:, 1] == ds]
    alpha = math.atan(CASES[case][4])
    c, s = math.cos(alpha), math.sin(alpha)
    # rows (t, id, x, y, z, qw, qx, qy, qz): the angle from the quaternion,
    # unwrapped as the reference angle
    theta = np.unwrap(2 * np.arctan2(q[:, 8], q[:, 5])) + alpha
    return (
        q[:, 0],
        q[:, 2] * c - q[:, 3] * s,
        q[:, 2] * s + q[:, 3] * c,
        theta,
        v[:, 2] * c - v[:, 3] * s,
        v[:, 2] * s + v[:, 3] * c,
        v[:, 7],
    )


def aos_command(program, case):
    """the command of the reference, its columns t x y theta vx vy w"""
    N, Rc, mu, mu_r, tan, vx0, vy0, w0, T = CASES[case]
    args = [N, Rc, mu, mu_r, math.atan(tan), vx0, vy0, w0, T, step(case), 0.5]
    return [program] + [str(a) if isinstance(a, int) else repr(float(a)) for a in args]


if len(sys.argv) == 4 and sys.argv[1] == "--vnative":
    run_vnative(sys.argv[2], sys.argv[3])
    sys.exit(0)

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import case_processes  # noqa: E402

program = os.path.abspath(sys.argv[1])
tmpdir = tempfile.mkdtemp()
files = {case: os.path.join(tmpdir, case + ".hdf5") for case in CASES}
# all the runs side by side: the vnative ones and the references
runs = case_processes.run(
    [[sys.executable, __file__, "--vnative", case, files[case]] for case in CASES]
    + [aos_command(program, case) for case in CASES]
)
vnative_runs = dict(zip(CASES, runs[: len(CASES)]))
refs = {}
failed = False
for case, r in zip(CASES, runs[len(CASES) :]):
    if r.returncode != 0:
        print(f"FAILED: {case}: reference exit status {r.returncode}\n{r.stderr}")
        failed = True
    else:
        refs[case] = np.loadtxt(r.stdout.splitlines())
for case in refs:
    ref = refs[case]
    r = vnative_runs[case]
    if r.returncode != 0:
        print(f"FAILED: {case}: exit status {r.returncode}\n{r.stdout}{r.stderr}")
        failed = True
        continue
    t, *traj = vnative(case, files[case])
    rows = np.rint(t / step(case)).astype(int)
    if rows.max() >= len(ref):
        print(f"FAILED: {case}: {len(t)} vnative rows, {len(ref)} reference rows")
        failed = True
        continue
    r = ref[rows]
    # x, y, theta then vx, vy, w; float32 tables: a few ulps of the largest
    # value. The contact of a vertex reaching the line can be activated one
    # step apart (the activation rules of the two kernels differ, with or
    # without rolling resistance): the velocities differ on a few rows at
    # each landing (the octagon with mu_r rocks: several landings per face),
    # and agree again after, the positions do not see it. So the positions
    # on every row, the velocities on all but 1% of the rows.
    errors = [np.abs(a - r[:, k + 1]) for k, a in enumerate(traj)]
    scale_q = max(1.0, np.max(np.abs(r[:, 1:4])))
    scale_v = max(1.0, np.max(np.abs(r[:, 4:7])))
    eq = np.max(errors[:3])
    ev = np.max(errors[3:], axis=0)
    off = np.count_nonzero(ev > 1e-5 * scale_v)
    ok = eq <= 1e-5 * scale_q and off <= 0.01 * len(t)
    x, y, theta = traj[0][-1], traj[1][-1], traj[2][-1]
    print(
        f"{case}: {len(t)} rows, max |q - q_ref| {eq:.2e}, max |v - v_ref| "
        f"{np.max(ev):.2e} ({off} rows above {1e-5 * scale_v:.0e}); (x, y, theta)(T) "
        f"({x:.6f}, {y:.6f}, {theta:.6f}) ref ({r[-1, 1]:.6f}, {r[-1, 2]:.6f}, "
        f"{r[-1, 3]:.6f}){'' if ok else ' FAILED'}"
    )
    failed |= not ok

# the physics, on the reference: the slide stops flat at v0^2 / (2 mu g),
# the square tips over without rolling resistance and falls back with it
ref = refs["slide"]
x_stop = 2.0**2 / (2 * 0.3 * g)
if (
    abs(ref[-1, 1] - x_stop) > 2 * 2.0 * step("slide")
    or np.max(np.abs(ref[:, 3])) > 1e-9
):
    print(f"FAILED: slide: x(T) {ref[-1, 1]:.6f}, expected {x_stop:.6f}, flat")
    failed = True
for case, over in (("tip-0", True), ("tip-0.1", False)):
    theta_T = refs[case][-1, 3]
    if (theta_T < -math.pi / 4) != over:
        print(f"FAILED: {case}: theta(T) {theta_T:.4f}, tips over expected {over}")
        failed = True

sys.exit(1 if failed else 0)
