#!/usr/bin/env python3
"""Validation of the TimeSteppingSchemesTest trajectories with hots-python.

hots-python is the stand-alone Python port of the Scilab code Hots (time-stepping schemes
for Lagrangian systems with frictionless unilateral contacts and a Newton impact law,
Acary, CMAME 2013). This script integrates, with hots-python, the test systems it can
model, on the time grid of the Siconos trajectories, and prints the maximal difference on
(q, v).

Usage:
    # 1. dump the full trajectories of the tests (in the build directory of the test)
    SICONOS_TEST_DUMP_TRAJECTORIES=1 ./TimeSteppingSchemesTest
    # 2. compare
    python3 validate_with_hots.py --hots <path to hots-python> [--dir <dir of the .dat>]

Conventions: hots-python is used with the Siconos conventions (see
docs/devel/PROJECTION_SCHEMES.md, section 5): activation of a contact for the
velocity-level problem when y_k + gamma h G(q_k) v_k <= 0 with gamma = 1/2, Jacobian G(q_k)
in the velocity-level problem. For linear constraints these are the hots-python rules;
for the nonlinear ones (rocking block), the Moreau step is rewritten below with the hots
building blocks.

Not covered (not in hots-python): friction, bilateral constraints, Newton-Euler systems,
first-order systems, nonlinear forces (the Moreau step of hots-python is linearized, the
Siconos one is solved by Newton iterations), D1Minus and MoreauJeanBilbao.
"""

import argparse
import math
import os
import sys

import numpy as np

parser = argparse.ArgumentParser(description=__doc__.split("\n")[0])
parser.add_argument(
    "--hots", default=os.environ.get("HOTS_PYTHON"), help="hots-python dir"
)
parser.add_argument(
    "--dir", default=".", help="directory of the <test>.dat trajectories"
)
parser.add_argument("--tol", type=float, default=1e-10, help="tolerance on (q, v)")
args = parser.parse_args()
if not args.hots:
    sys.exit("give the path of hots-python (--hots or HOTS_PYTHON)")
sys.path.insert(0, args.hots)

from hots.model import Model, Trajectory
from hots.moreau import _free_velocity, _impulse_step
from hots.paoli_schatzman import paoli_schatzman_fixed
from hots.projected import (
    combined_projected_moreau_fixed,
    projected_moreau_fixed,
)

g = 9.81


def zero_kc(n):
    return lambda t, q, v: (np.zeros((n, n)), np.zeros((n, n)))


def zero_fint(n):
    return lambda t, q, v: np.zeros(n)


def moreau_siconos(t, model, theta=0.5, gamma=0.5):
    """Moreau-Jean with the Siconos conventions: activation y_k + gamma h G(q_k) v_k <= 0
    (linear forecast) and G(q_k) in the velocity-level problem."""
    n = t.size
    q = np.zeros((n, model.ndof))
    v = np.zeros((n, model.ndof))
    q[0], v[0] = model.q0, model.v0
    status = model.status0
    for k in range(n - 1):
        h = t[k + 1] - t[k]
        W, vfree = _free_velocity(model, t[k], t[k + 1], q[k], v[k], theta)
        G = model.G(q[k])
        active = np.flatnonzero(model.H(q[k]) + gamma * h * G @ v[k] <= 0.0)
        s = _impulse_step(
            model, t[k], t[k + 1], q[k], v[k], status, W, vfree, active, G, theta, 1e-10
        )
        q[k + 1], v[k + 1], status = s.q, s.v, s.status
    return Trajectory(t, q, v)


def paoli_with_kc(t, model, theta=0.5):
    """Schatzman-Paoli scheme of SchatzmanPaoliOSI for M q'' + C q' + K q = F (constant M,
    K, C, F). With D2q = q_{k+1} - 2 q_k + q_{k-1}:

        M D2q + h^2 C (theta (q_{k+1} - q_{k-1}) / (2h) + (1 - theta) v_k)
              + h^2 K (theta^2 q_{k+1} + (1 - theta^2) q_k) - h^2 F = G^T r,
        0 <= g(q_{k+1}) + e g(q_{k-1})  _|_  r >= 0,   v_{k+1} = (q_{k+1} - q_{k-1}) / (2h),

    started from q_{-1} = q_0 - h v_0 + h^2/2 M^{-1} (F - K q_0 - C v_0).
    hots.paoli_schatzman (a port of the Scilab PaoliLinearOnestep) uses other weights of K
    and C and is much less accurate on this system (3e-3 instead of 7e-5 before the first
    impact of springChainSchatzmanPaoli, against the exact solution), so it is not used for
    systems with K or C."""
    from hots.lcp import lcp

    n = t.size
    h = t[1] - t[0]
    K, C = model.KC(t[0], model.q0, model.v0)
    M = model.M(model.q0, t[0])
    F = model.fext(t[0])
    Winv = np.linalg.inv(M + 0.5 * h * theta * C + h * h * theta * theta * K)
    qm = (
        model.q0
        - h * model.v0
        + 0.5 * h * h * np.linalg.solve(M, F - K @ model.q0 - C @ model.v0)
    )
    qi, vi = model.q0.copy(), model.v0.copy()
    q = np.zeros((n, model.ndof))
    v = np.zeros((n, model.ndof))
    q[0], v[0] = qi, vi
    for k in range(n - 1):
        residu = (
            M @ (qm - qi)
            + h * h * C @ (theta * (qi - qm) / (2 * h) + (1 - theta) * vi)
            + h * h * K @ qi
            - h * h * F
        )
        qfree = qi - Winv @ residu
        G = model.G(qi)
        r, _ = lcp(G @ Winv @ G.T, model.H(qfree) + model.e * model.H(qm))
        qn = qfree + Winv @ G.T @ r
        vi = (qn - qm) / (2 * h)
        qm, qi = qi, qn
        q[k + 1], v[k + 1] = qi, vi
    return Trajectory(t, q, v)


def combined_as_siconos(t, model, theta=0.5, tolactive=1e-10):
    """Combined projection (Algorithm 2) as computed by TimeSteppingCombinedProjection: the
    constraints to add to the active set are those violated at a = q_k + h (theta v_{k+1}
    + (1 - theta) v_k), before the projection on the active set (Hots/CMAME: after it)."""
    from hots.lcp import lcp

    n = t.size
    q = np.zeros((n, model.ndof))
    v = np.zeros((n, model.ndof))
    q[0], v[0] = model.q0, model.v0
    for k in range(n - 1):
        h = t[k + 1] - t[k]
        qi, vi = q[k], v[k]
        W, vfree = _free_velocity(model, t[k], t[k + 1], qi, vi, theta)
        Gi = model.G(qi)
        active = []
        while True:
            vf = vfree
            if active:
                Ga = Gi[active]
                r, _ = lcp(Ga @ W @ Ga.T, Ga @ vfree + model.e * Ga @ vi)
                vf = vfree + W @ Ga.T @ r
            a = qi + h * (theta * vf + (1 - theta) * vi)
            qf = a
            if active:
                for _ in range(101):
                    qold = qf
                    Gf = model.G(qold)[active]
                    bf = model.H(qold)[active] + Gf @ (a - qold)
                    qf = a + Gf.T @ np.linalg.solve(Gf @ Gf.T, -bf)
                    if np.linalg.norm(qf - qold) <= 1e-12:
                        break
            violated = [
                int(j)
                for j in np.flatnonzero(model.H(a) < -tolactive)
                if j not in active
            ]
            if not violated:
                break
            active.extend(violated)
        q[k + 1], v[k + 1] = qf, vf
    return Trajectory(t, q, v)


# ---------------------------------------------------------------------------------------
# Models (data of TimeSteppingSchemesTest.cpp)
# ---------------------------------------------------------------------------------------
def beads():
    R = 0.1

    def h(q):
        return np.array([q[0] - R, q[1] - q[0] - 2 * R, q[2] - q[1] - 2 * R])

    G = np.array([[1.0, 0, 0], [-1, 1, 0], [0, -1, 1]])
    return Model(
        name="beads",
        q0=[0.15, 0.45, 0.8],
        v0=[0.0, -1.0, -2.0],
        mass=lambda q, t: np.eye(3),
        fint=zero_fint(3),
        fext=lambda t: np.full(3, -g),
        nabla_fint=zero_kc(3),
        h=h,
        nabla_h=lambda q: G,
        e=0.8,
    )


def rocking_block():
    L, Hb, m = 0.2, 0.1, 1.0
    th0 = 0.3
    q0 = [
        0.5 * L * math.cos(th0) - 0.5 * Hb * math.sin(th0),
        0.5 * L * math.sin(th0) + 0.5 * Hb * math.cos(th0),
        th0,
    ]
    M = np.diag([m, m, m * (L * L + Hb * Hb) / 12])

    def h(q):
        c, s = math.cos(q[2]), math.sin(q[2])
        return np.array(
            [q[1] - 0.5 * L * s - 0.5 * Hb * c, q[1] + 0.5 * L * s - 0.5 * Hb * c]
        )

    def G(q):
        c, s = math.cos(q[2]), math.sin(q[2])
        return np.array(
            [[0, 1, -0.5 * L * c + 0.5 * Hb * s], [0, 1, 0.5 * L * c + 0.5 * Hb * s]]
        )

    return Model(
        name="rocking_block",
        q0=q0,
        v0=np.zeros(3),
        mass=lambda q, t: M,
        fint=zero_fint(3),
        fext=lambda t: np.array([0.0, -m * g, 0.0]),
        nabla_fint=zero_kc(3),
        h=h,
        nabla_h=G,
        e=0.5,
        is_linear=False,
    )


def spring_chain(clamped):
    """Chain of 5 masses; if clamped, the first one is removed (q_0 = v_0 = 0)."""
    n, k, damping = 5, 100.0, 0.01
    K = np.diag([2 * k] * 4 + [k]) - np.diag([k] * 4, 1) - np.diag([k] * 4, -1)
    dofs = list(range(1, n)) if clamped else list(range(n))
    K = K[np.ix_(dofs, dofs)]
    C = damping * K
    nd = len(dofs)
    i4, i2 = dofs.index(4), dofs.index(2)

    def h(q):
        return np.array([0.05 - q[i4], 0.08 - q[i2]])

    G = np.zeros((2, nd))
    G[0, i4], G[1, i2] = -1.0, -1.0
    return Model(
        name="chain",
        q0=np.zeros(nd),
        v0=np.ones(nd),
        mass=lambda q, t: np.eye(nd),
        fint=lambda t, q, v: K @ q + C @ v,
        fext=lambda t: np.zeros(nd),
        nabla_fint=lambda t, q, v: (K, C),
        h=h,
        nabla_h=lambda q: G,
        e=0.5,
    )


# ---------------------------------------------------------------------------------------
# Siconos columns of (q, v): [t, q_ds, v_ds, ...] for each dynamical system
# ---------------------------------------------------------------------------------------
def beads_qv(d):
    return d[:, [1, 3, 5]], d[:, [2, 4, 6]]


def block_qv(d):
    return d[:, 1:4], d[:, 4:7]


def chain_qv(clamped):
    def f(d):
        q, v = d[:, 1:6], d[:, 6:11]
        return (q[:, 1:], v[:, 1:]) if clamped else (q, v)

    return f


def run_hots(scheme, model, t):
    if scheme == "moreau":
        return moreau_siconos(t, model)
    if scheme == "combined":
        # Siconos: setConstraintTolUnilateral(1e-10), exact projection (from_origin)
        return combined_projected_moreau_fixed(
            t, model, theta=0.5, tolactive=1e-10, fpitermax=100, fptol=1e-12
        )
    if scheme == "projected":
        return projected_moreau_fixed(
            t, model, theta=0.5, gamma=0.5, tolactive=1e-7, fpitermax=100, fptol=1e-12
        )
    if scheme == "paoli":
        return paoli_schatzman_fixed(t, model, theta=0.5)
    if scheme == "paoli-kc":
        return paoli_with_kc(t, model, theta=0.5)
    if scheme == "combined-as-siconos":
        return combined_as_siconos(t, model)
    raise ValueError(scheme)


# test name: (model, columns of (q, v), hots scheme, compare v, tolerance)
# A name may be listed twice (trailing space) to compare it with a second hots scheme.
cases = {
    "beadsMoreauJean": (beads, beads_qv, "moreau", True, None),
    "beadsManagedMoreauJean": (beads, beads_qv, "moreau", True, None),
    "beadsMoreauJeanSparseBlockStorage": (beads, beads_qv, "moreau", True, None),
    "beadsDirectProjection": (beads, beads_qv, "projected", True, None),
    # Known deviation: Siconos tests the violation before the projection (see
    # combined_as_siconos), Hots/CMAME after it.
    "beadsCombinedProjection": (beads, beads_qv, "combined", True, None),
    "beadsCombinedProjection ": (beads, beads_qv, "combined-as-siconos", True, None),
    # v is a difference quotient in this scheme, it is not compared.
    "beadsSchatzmanPaoli": (beads, beads_qv, "paoli", False, None),
    "rockingBlockMoreauJean": (rocking_block, block_qv, "moreau", True, None),
    # the projections stop when the gaps are below 1e-10 (Siconos) or when the fixed point
    # converges (Hots, 1e-12): the positions differ by a few 1e-10
    "rockingBlockCombinedProjectionFromOrigin": (
        rocking_block,
        block_qv,
        "combined",
        True,
        1e-9,
    ),
    "springChainClampedMoreauJean": (
        lambda: spring_chain(True),
        chain_qv(True),
        "moreau",
        True,
        None,
    ),
    "springChainMoreauJean": (
        lambda: spring_chain(False),
        chain_qv(False),
        "moreau",
        True,
        None,
    ),
    # Siconos formulas transcribed in Python (see paoli_with_kc); v is a difference quotient.
    "springChainSchatzmanPaoli": (
        lambda: spring_chain(False),
        chain_qv(False),
        "paoli-kc",
        False,
        None,
    ),
}

failed = 0
for name, (model_fn, qv, scheme, with_v, case_tol) in cases.items():
    tol = case_tol or args.tol
    path = os.path.join(args.dir, name.strip() + ".dat")
    if not os.path.exists(path):
        print(f"{name:45s} no trajectory {path}")
        continue
    d = np.loadtxt(path, ndmin=2)
    t = np.concatenate([[0.0], d[:, 0]])
    tr = run_hots(scheme, model_fn(), t)
    q, v = qv(d)
    err_q = np.abs(q - tr.q[1:]).max()
    err_v = np.abs(v - tr.v[1:]).max() if with_v else 0.0
    err = max(err_q, err_v)
    ok = err <= tol
    failed += not ok
    k = int(np.argmax(np.abs(q - tr.q[1:]).max(axis=1) > tol)) if not ok else -1
    print(
        f"{name.strip():41s} {scheme:20s} max |q - q_hots| {err_q:9.2e}"
        + (f"  max |v - v_hots| {err_v:9.2e}" if with_v else " " * 28)
        + ("  ok" if ok else f"  DIFFERENT (from t = {d[k, 0]:.3f})")
    )
sys.exit(1 if failed else 0)
