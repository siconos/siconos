"""Python port of the vkernel 2D contact manifold (cxd_narrowphase.hpp).

normal  : separated -> closest-feature direction A -> B (distance GJK)
          penetrating -> least-penetration axis (EPA == SAT for polygons)
manifold: clip_manifold (reference / incident face clipping, face_tol 1e-6,
          reference hysteresis 1e-3, ties to A), points sorted along the
          tangent (-n_y, n_x); no face contact -> the single closest point.
"""

import numpy as np

FACE_TOL = 1e-6
INACTIVE_GAP = 1e30


def _edges(V):
    n = len(V)
    return [(i, (i + 1) % n) for i in range(n)] if n >= 3 else [(0, 1)]


def sat(A, B):
    best, best_n = -np.inf, None
    for V in (A, B):
        for i, j in _edges(V):
            e = V[j] - V[i]
            u = np.array([e[1], -e[0]]) / np.linalg.norm(e)
            pa, pb = A @ u, B @ u
            s_pos, s_neg = pb.min() - pa.max(), pa.min() - pb.max()
            sep = max(s_pos, s_neg)
            if sep > best:
                best, best_n = sep, (u if s_pos >= s_neg else -u)
    return best, best_n


def _foot(p, a, b):
    e = b - a
    t = np.clip(np.dot(p - a, e) / np.dot(e, e), 0, 1)
    return a + t * e


def closest(A, B):
    best = (np.inf, None, None)
    for P, Q, flip in ((A, B, False), (B, A, True)):
        for p in P:
            for i, j in _edges(Q):
                f = _foot(p, Q[i], Q[j])
                d = np.linalg.norm(f - p)
                if d < best[0]:
                    best = (d, f, p) if flip else (d, p, f)
    return best  # distance, point on A, point on B


def _best_face(V, d):
    nv = len(V)
    area2 = sum(
        V[i][0] * V[(i + 1) % nv][1] - V[(i + 1) % nv][0] * V[i][1] for i in range(nv)
    )
    orient = 1.0 if area2 >= 0 else -1.0
    best = None
    for i in range(nv):
        j = (i + 1) % nv
        e = V[j] - V[i]
        o = np.array([e[1], -e[0]]) / np.linalg.norm(e)
        if nv >= 3:
            o = orient * o
        elif np.dot(o, d) < 0:
            o = -o
        al = np.dot(o, d)
        if best is None or al > best[3]:
            best = (i, j, o, al)
    return best


def clip(A, B, n):
    fa, fb = _best_face(A, n), _best_face(B, -n)
    if max(fa[3], fb[3]) < 1 - FACE_TOL:
        return None
    ref_a = not (fb[3] > fa[3] + 1e-3)
    R, I = (A, B) if ref_a else (B, A)
    rf, inc = (fa, fb) if ref_a else (fb, fa)
    r0, r1 = R[rf[0]], R[rf[1]]
    p0, p1 = I[inc[0]], I[inc[1]]
    L = np.linalg.norm(r1 - r0)
    t = (r1 - r0) / L
    s0, s1 = np.dot(t, p0 - r0), np.dot(t, p1 - r0)
    if abs(s1 - s0) <= 1e-12 * (1 + L):
        if s0 < 0 or s0 > L:
            return None
        lmin, lmax = 0.0, 1.0
    else:
        la, lb = -s0 / (s1 - s0), (L - s0) / (s1 - s0)
        lmin, lmax = max(0.0, min(la, lb)), min(1.0, max(la, lb))
        if lmin > lmax:
            return None
    dp = p1 - p0
    q = [p0 + lmin * dp, p0 + lmax * dp]
    nq = 2 if (lmax - lmin) ** 2 * np.dot(dp, dp) > (1e-12 * (1 + L)) ** 2 else 1
    pts = []
    for k in range(nq):
        gap = np.dot(rf[2], q[k] - r0)
        w = q[k] - gap * rf[2]
        pts.append(dict(h=gap, cp1=w if ref_a else q[k], cp2=q[k] if ref_a else w))
    if nq == 2:
        tng = np.array([-n[1], n[0]])
        if np.dot(tng, pts[1]["cp1"]) < np.dot(tng, pts[0]["cp1"]):
            pts.reverse()
    info = dict(ref_a=ref_a, ref=(r0, r1), inc=(p0, p1), ref_normal=rf[2])
    return pts, info


def manifold(A, B):
    """A, B: world vertices (polygon CCW, or a 2-vertex segment)."""
    A, B = np.asarray(A, float), np.asarray(B, float)
    sep, n_sat = sat(A, B)
    if sep > 0:
        d, pa, pb = closest(A, B)
        n, h, cp1, cp2 = (pb - pa) / d, d, pa, pb
    else:
        n, h = n_sat, sep
        cp1 = cp2 = None
    c = clip(A, B, n)
    if c is None:
        if cp1 is None:  # penetrating vertex-vertex: deepest vertex witnesses
            cp1 = A[np.argmax(A @ n)]
            cp2 = B[np.argmin(B @ n)]
        return dict(normal=n, h=h, points=[dict(h=h, cp1=cp1, cp2=cp2)], info=None)
    pts, info = c
    return dict(normal=n, h=h, points=pts, info=info)


def point_gaps(m):
    """gap of manifold point 0 and 1 (a missing point: INACTIVE_GAP)."""
    g = [p["h"] for p in m["points"]]
    return g + [INACTIVE_GAP] * (2 - len(g))


def rect(w, h):
    return np.array(
        [[-w / 2, -h / 2], [w / 2, -h / 2], [w / 2, h / 2], [-w / 2, h / 2]]
    )


def place(V, x, y, th):
    c, s = np.cos(th), np.sin(th)
    return V @ np.array([[c, s], [-s, c]]) + np.array([x, y])
