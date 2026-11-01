"""Read a scene's hdf5 output and rebuild, per step, the manifold of every
polytope pair with manifold.py, matched to the simulator's active rows."""

import h5py
import numpy as np

from manifold import manifold, place, point_gaps, rect

H = 1e-3
SCENES = {
    "flat": dict(bodies={1: rect(1.0, 0.4)}, segs=[(-3, 0, 3, 0)]),
    "tilted": dict(bodies={1: rect(1.0, 0.4)}, segs=[(-3, 0, 3, 0)]),
    "overhang": dict(bodies={1: rect(1.0, 0.4)}, segs=[(-3, 0, 0, 0)]),
    "tipping": dict(bodies={1: rect(1.0, 0.4)}, segs=[(-3, 0, 0, 0)]),
    "junction": dict(bodies={1: rect(1.0, 0.4)}, segs=[(-3, 0, 0, 0), (0, 0, 3, 0)]),
    "stack": dict(bodies={1: rect(2.0, 0.4), 2: rect(0.6, 0.3)}, segs=[(-3, 0, 3, 0)]),
    "corner": dict(bodies={1: rect(2.0, 0.4), 2: rect(0.4, 0.4)}, segs=[(-3, 0, 3, 0)]),
}


def load(scene):
    with h5py.File(f"{scene}.hdf5") as f:
        dyn = f["data"]["dynamic"][:]
        cf = f["data"]["cf"][:]
        vel = f["data"]["velocities"][:]
    ts = np.unique(dyn[:, 0])
    q, v = {}, {}
    for bid in SCENES[scene]["bodies"]:
        d = dyn[dyn[:, 1] == bid]
        q[bid] = np.c_[d[:, 2], d[:, 3], 2 * np.arctan2(d[:, 8], d[:, 5])]
        w = vel[vel[:, 1] == bid]
        v[bid] = np.c_[w[:, 2], w[:, 3], w[:, 7]]
    return ts, q, v, cf


def world(scene, q, bid, k):
    x, y, th = q[bid][k]
    return place(SCENES[scene]["bodies"][bid], x, y, th)


def pairs(scene):
    """(label, ds1 body id, other: body id or segment tuple)"""
    sc = SCENES[scene]
    out = []
    ids = sorted(sc["bodies"])
    for b in ids:
        for s in sc["segs"]:
            out.append((f"body{b}-seg{s}", b, s))
    for i, a in enumerate(ids):
        for b in ids[i + 1 :]:
            out.append((f"body{a}-body{b}", a, b))
    return out


def manifolds(scene, q, k):
    res = {}
    for label, a, other in pairs(scene):
        A = world(scene, q, a, k)
        B = (
            np.array([[other[0], other[1]], [other[2], other[3]]], float)
            if isinstance(other, tuple)
            else world(scene, q, other, k)
        )
        res[label] = manifold(A, B)
    return res


def cf_rows(cf, t):
    """active rows at time t: dict(cp1, cp2, n (1->2), y, lam, ids)"""
    rows = cf[np.isclose(cf[:, 0], t)]
    return [
        dict(
            cp2=r[2:4],
            cp1=r[5:7],
            n=-r[8:10],
            y=r[14],
            ydot=r[17],
            lam_n=r[20],
            lam_t=r[21],
            id1=int(r[24]),
            id2=int(r[25]),
        )
        for r in rows
    ]


def match(scene, mans, row):
    """the (pair label, point index) of an active row: nearest cp1"""
    best = (np.inf, None, None)
    for label, m in mans.items():
        _, a, other = next(p for p in pairs(scene) if p[0] == label)
        if a != row["id1"] and a != row["id2"]:
            continue
        for i, p in enumerate(m["points"]):
            d = min(
                np.linalg.norm(p["cp1"] - row["cp1"])
                + np.linalg.norm(p["cp2"] - row["cp2"]),
                np.linalg.norm(p["cp2"] - row["cp1"])
                + np.linalg.norm(p["cp1"] - row["cp2"]),
            )
            if d < best[0]:
                best = (d, label, i)
    return best


def point_velocity(v, q, bid, k, p):
    vx, vy, w = v[bid][k]
    r = p - q[bid][k, :2]
    return np.array([vx - w * r[1], vy + w * r[0]])


def rate(scene, q, v, k, label, pt, n):
    """d/dt of the gap n.(cp2 - cp1) at configuration k"""
    _, a, other = next(p for p in pairs(scene) if p[0] == label)
    va = point_velocity(v, q, a, k, pt["cp1"])
    vb = (
        np.zeros(2)
        if isinstance(other, tuple)
        else point_velocity(v, q, other, k, pt["cp2"])
    )
    return float(np.dot(n, vb - va))


def history(scene, gamma=0.5):
    """per pair and manifold point: gap at the end of each step, the
    activation predicted from the start of the step (y + gamma h ydot <= 0),
    the simulator's activation (cf rows) and normal impulse"""
    ts, q, v, cf = load(scene)
    labels = [p[0] for p in pairs(scene)]
    z = lambda: np.zeros((len(ts), 2))
    hist = {
        l: dict(
            gap=np.full((len(ts), 2), np.nan),
            pred=z().astype(bool),
            active=z().astype(bool),
            lam=z(),
            lamt=z(),
        )
        for l in labels
    }
    mans_k = [manifolds(scene, q, k) for k in range(len(ts))]
    worst, mismatch = 0.0, 0
    for k, t in enumerate(ts):
        for l, m in mans_k[k].items():
            hist[l]["gap"][k] = point_gaps(m)
        if k == 0:
            continue
        # the step (k-1 -> k) is decided at configuration k-1
        cand = []
        for l, m in mans_k[k - 1].items():
            for i, pt in enumerate(m["points"]):
                yd = rate(scene, q, v, k - 1, l, pt, m["normal"])
                if pt["h"] + gamma * H * yd <= 1e-7:  # hdf5 stores float32
                    hist[l]["pred"][k, i] = True
                    cand.append((l, i, pt))
        rows = cf_rows(cf, t)
        if len(rows) != len(cand):
            mismatch += 1
        # exclusive nearest assignment of the rows to the predicted points
        free = list(cand)
        for row in rows:
            if not free:
                break
            ds = [
                min(
                    np.linalg.norm(pt["cp1"] - row["cp1"])
                    + np.linalg.norm(pt["cp2"] - row["cp2"]),
                    np.linalg.norm(pt["cp2"] - row["cp1"])
                    + np.linalg.norm(pt["cp1"] - row["cp2"]),
                )
                for (_, _, pt) in free
            ]
            j = int(np.argmin(ds))
            l, i, _ = free.pop(j)
            worst = max(worst, ds[j])
            hist[l]["active"][k, i] = True
            hist[l]["lam"][k, i] = row["lam_n"]
            hist[l]["lamt"][k, i] = row["lam_t"]
    return ts, q, v, cf, hist, worst, mismatch
