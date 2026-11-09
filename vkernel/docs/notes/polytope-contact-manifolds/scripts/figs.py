"""Figures of the polytope contact manifolds (vkernel, Moreau-Jean)."""

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.patches import Polygon

from manifold import INACTIVE_GAP, manifold, place, rect
from sim import SCENES, history, manifolds, pairs

C0, C1, C2, C3 = "#2a78d6", "#eb6834", "#1baf7a", "#eda100"  # point 0 / 1, 2nd pair
INK, MUTED, BODY = "#0b0b0b", "#52514e", "#d9d8d4"
plt.rcParams.update(
    {
        "font.size": 9,
        "axes.edgecolor": MUTED,
        "axes.labelcolor": INK,
        "xtick.color": MUTED,
        "ytick.color": MUTED,
        "axes.grid": True,
        "grid.color": "#e6e5e1",
        "grid.linewidth": 0.6,
        "figure.facecolor": "#fcfcfb",
        "axes.facecolor": "#fcfcfb",
    }
)
MARK = ["o", "s"]


def draw_body(ax, V):
    ax.add_patch(Polygon(V, closed=True, fc=BODY, ec=MUTED, lw=1.0, zorder=1))


def draw_seg(ax, s, color=INK):
    ax.plot(
        [s[0], s[2]], [s[1], s[3]], color=color, lw=2.5, solid_capstyle="butt", zorder=2
    )
    ax.plot([s[0], s[2]], [s[1], s[3]], "|", color=color, ms=8, mew=1.5, zorder=2)


def draw_points(
    ax,
    m,
    colors=(C0, C1),
    active=(None, None),
    lam=(0, 0),
    tag="",
    annotate=True,
    lift=0,
):
    """manifold points: cp on the incident shape (marker) and its witness on
    the reference face (x), joined by the gap; filled = active"""
    for i in range(2):
        if i >= len(m["points"]):
            if annotate:
                ax.text(
                    0.98,
                    0.04 + 0.08 * (1 - i),
                    f"{tag}p{i}: no point (h = 1e30)",
                    transform=ax.transAxes,
                    ha="right",
                    color=MUTED,
                    fontsize=8,
                )
            continue
        p = m["points"][i]
        a, b = p["cp1"], p["cp2"]
        ax.plot([a[0], b[0]], [a[1], b[1]], color=colors[i], lw=1, zorder=3)
        is_on = active[i] if active[i] is not None else p["h"] <= 0
        size = 6 + (min(14, 6 * np.sqrt(lam[i] / 200)) if lam[i] > 0 else 0)
        for c in (a, b):
            ax.plot(
                *c,
                MARK[i],
                ms=size,
                mfc=colors[i] if is_on else "#fcfcfb",
                mec=colors[i],
                mew=1.6,
                zorder=4,
            )
        if annotate:
            h = p["h"]
            txt = f"{tag}p{i}  h={h:+.3f}" if abs(h) >= 5e-4 else f"{tag}p{i}  h≈0"
            ax.annotate(
                txt,
                b,
                xytext=(4, (-12 if i == 0 else 6) + lift),
                textcoords="offset points",
                color=INK,
                fontsize=8,
                zorder=5,
            )


def draw_construction(ax, m):
    """reference face (thick), incident face (dashed) and the clipping lines"""
    if m["info"] is None:
        return
    r0, r1 = m["info"]["ref"]
    p0, p1 = m["info"]["inc"]
    nr = m["info"]["ref_normal"]
    ax.plot(
        [r0[0], r1[0]],
        [r0[1], r1[1]],
        color=C2,
        lw=5,
        alpha=0.45,
        zorder=2,
        solid_capstyle="butt",
        label="reference face",
    )
    ax.plot(
        [p0[0], p1[0]],
        [p0[1], p1[1]],
        color=C3,
        lw=2,
        ls="--",
        zorder=2,
        label="incident face",
    )
    for r in (r0, r1):
        e = np.array([r - 0.6 * nr, r + 0.6 * nr])
        ax.plot(e[:, 0], e[:, 1], ":", color=MUTED, lw=1, zorder=2)
    c = 0.5 * (r0 + r1)
    ax.annotate(
        "",
        c + 0.18 * nr,
        c,
        arrowprops=dict(arrowstyle="->", color=C2, lw=1.5),
        zorder=5,
    )


def vertex_at(V, p):
    """V translated so that its lowest vertex is at p"""
    return V - V[np.argmin(V[:, 1])] + np.asarray(p)


def fig_cases(path):
    seg = (-1.0, 0.0, 1.0, 0.0)
    S = np.array([[seg[0], seg[1]], [seg[2], seg[3]]])
    B = rect(1.0, 0.4)
    th = np.radians(15)
    cases = [
        (
            "(a) face on face, flat",
            [place(B, 0.0, 0.2 + 0.03, 0)],
            [S],
            "parallel faces: tie -> body 1 (box) is the\nreference; the segment is clipped to it",
        ),
        (
            "(b) corner first, tilted 15°",
            [place(B, 0.0, 0.2 * np.cos(th) + 0.5 * np.sin(th) + 0.01, th)],
            [S],
            "both corners of the lower face: p1 is\nthe far corner, gap L·sin θ (inactive)",
        ),
        (
            "(c) overhang: box past the segment end",
            [place(B, 0.75, 0.2 + 0.02, 0)],
            [S],
            "p1 is clipped to the segment END\n(not a box corner)",
        ),
        (
            "(d) tipping over the segment end",
            [place(B, 1.25, 0.11, -0.45)],
            [S],
            "vertex (segment end) on the box face:\nreference = box face",
        ),
        (
            "(e) two collinear segments (mesh edges)",
            [place(B, 0.1, 0.2 + 0.02, 0)],
            [np.array([[-1.0, 0], [0, 0]]), np.array([[0, 0], [1.0, 0]])],
            "one manifold per segment: 4 points,\ntwo of them at the junction",
        ),
        (
            "(f) vertex on vertex",
            [vertex_at(place(B, 0, 0, np.radians(45)), (1.03, 0.12))],
            [S],
            "no face aligned with the normal\n(1e-6): ONE point, p1 gets h = 1e30",
        ),
    ]
    fig, axs = plt.subplots(2, 3, figsize=(13, 7.2))
    for ax, (title, bodies, segs, note) in zip(axs.flat, cases):
        for V in bodies:
            draw_body(ax, V)
        for k, s in enumerate(segs):
            draw_seg(ax, (s[0][0], s[0][1], s[1][0], s[1][1]))
            m = manifold(bodies[0], s)
            draw_construction(ax, m)
            draw_points(
                ax,
                m,
                colors=(C0, C1) if k == 0 else (C2, C3),
                tag=f"s{k + 1}:" if len(segs) > 1 else "",
                lift=-16 * k if k else 0,
            )
            n = m["normal"]
            c = np.mean([p["cp1"] for p in m["points"]], axis=0)
            ax.annotate(
                "",
                c + 0.25 * n,
                c,
                arrowprops=dict(arrowstyle="-|>", color=INK, lw=1.2),
            )
        ax.set_title(title, loc="left", color=INK)
        ax.text(
            0.02, 0.97, note, transform=ax.transAxes, va="top", color=MUTED, fontsize=8
        )
        ax.set_aspect("equal")
        ax.set_xlim(
            -1.1, 1.9 if any(x in title for x in ("(c)", "(d)", "(f)")) else 1.1
        )
        ax.set_ylim(-0.35, 0.85)
    h, l = axs[0, 1].get_legend_handles_labels()
    fig.legend(h[:2], l[:2], loc="lower center", ncol=2, frameon=False)
    fig.suptitle(
        "Contact manifold of a polytope against a segment: reference / incident face clipping\n"
        "black arrow: contact normal (body 1 → body 2); dotted: clipping planes; "
        "filled marker: h ≤ 0, hollow: h > 0",
        color=INK,
    )
    fig.tight_layout(rect=(0, 0.04, 1, 1))
    fig.savefig(path, dpi=130)
    plt.close(fig)


def snapshot(ax, scene, q, hist, k, xlim, ylim, title):
    """configuration k-1 as the solver saw it for the step k-1 -> k"""
    kc = max(k - 1, 0)
    sc = SCENES[scene]
    for bid in sc["bodies"]:
        x, y, th = q[bid][kc]
        draw_body(ax, place(sc["bodies"][bid], x, y, th))
    for s in sc["segs"]:
        draw_seg(ax, s)
    mans = manifolds(scene, q, kc)
    cols = [(C0, C1), (C2, C3)]
    j = 0
    for label, m in mans.items():
        gaps = [p["h"] for p in m["points"]]
        if min(gaps) > 0.3:
            continue
        hh = hist[label]
        multi = len(sc["segs"]) > 1 or len(sc["bodies"]) > 1
        draw_points(
            ax,
            m,
            colors=pair_colors(label),
            active=tuple(hh["active"][k]),
            lam=tuple(hh["lam"][k]),
            tag=(short(label) + " ") if multi else "",
            annotate=True,
        )
        j += 1
    ax.set_aspect("equal")
    ax.set_xlim(*xlim)
    ax.set_ylim(*ylim)
    ax.set_title(title, loc="left", color=INK, fontsize=9)


def pair_colors(label):
    return (C2, C3) if short(label) in ("s2", "bb") else (C0, C1)


def short(label):
    return {
        "body1-seg(-3, 0, 0, 0)": "s1",
        "body1-seg(0, 0, 3, 0)": "s2",
        "body1-seg(-3, 0, 3, 0)": "g",
        "body1-body2": "bb",
        "body2-seg(-3, 0, 3, 0)": "g2",
    }.get(label, label)


def series(axg, axl, ts, hh, colors=(C0, C1), tag=""):
    for i in range(2):
        g = np.where(hh["gap"][:, i] >= INACTIVE_GAP / 2, np.nan, hh["gap"][:, i])
        axg.plot(ts, g, color=colors[i], lw=1.5, label=f"{tag}p{i} gap h")
        a = hh["active"][:, i]
        axg.plot(ts[a], np.full(a.sum(), -2e-4 - 4e-4 * i), "|", color=colors[i], ms=6)
        axl.plot(ts, hh["lam"][:, i], color=colors[i], lw=1.2, label=f"{tag}p{i}")
    axg.set_yscale("symlog", linthresh=1e-3)
    axg.axhline(0, color=MUTED, lw=0.8)
    axg.set_ylabel("gap h of the point [m]")
    axl.set_ylabel("normal impulse λn [N·s]")
    axl.set_xlabel("t [s]")


def fig_scene_ts(path, scene, snaps, xlim, ylim, label, tmax=None, title=""):
    ts, q, v, cf, hist, worst, mism = history(scene)
    fig = plt.figure(figsize=(13, 7.5))
    gs = fig.add_gridspec(3, len(snaps), height_ratios=[1.25, 1, 0.8])
    for c, t in enumerate(snaps):
        k = int(np.argmin(np.abs(ts - t)))
        snapshot(
            fig.add_subplot(gs[0, c]),
            scene,
            q,
            hist,
            k,
            xlim,
            ylim,
            f"t = {ts[k - 1]:.3f} s",
        )
    axg = fig.add_subplot(gs[1, :])
    axl = fig.add_subplot(gs[2, :], sharex=axg)
    series(axg, axl, ts, hist[label], colors=pair_colors(label), tag=short(label) + " ")
    for t in snaps:
        axg.axvline(t, color=MUTED, lw=0.6, ls=":")
    if tmax:
        axg.set_xlim(0, tmax)
    axg.legend(loc="upper right", frameon=False)
    axg.text(
        0.005,
        0.04,
        "ticks below 0: steps where the interaction is in the active set (simulator)",
        transform=axg.transAxes,
        color=MUTED,
        fontsize=8,
    )
    fig.suptitle(
        title
        + f"\n(active set predicted from the manifold = simulator's on {len(ts) - 1 - mism}/{len(ts) - 1} steps)",
        color=INK,
    )
    fig.tight_layout()
    fig.savefig(path, dpi=130)
    plt.close(fig)
    return hist, ts


def fig_static(path, rows):
    """one snapshot per scene, with the share of the normal impulse"""
    fig, axs = plt.subplots(1, len(rows), figsize=(4.4 * len(rows), 3.8))
    for ax, (scene, t, xlim, ylim, title) in zip(np.atleast_1d(axs), rows):
        ts, q, v, cf, hist, worst, mism = history(scene)
        k = int(np.argmin(np.abs(ts - t)))
        snapshot(ax, scene, q, hist, k, xlim, ylim, title)
        lines = []
        for label, hh in hist.items():
            for i in range(2):
                if hh["active"][k, i]:
                    lines.append(f"{short(label)} p{i}: λn={hh['lam'][k, i]:.2f}")
        ax.text(
            0.02,
            0.97,
            "\n".join(lines),
            transform=ax.transAxes,
            va="top",
            fontsize=8,
            color=INK,
        )
    fig.tight_layout()
    fig.savefig(path, dpi=130)
    plt.close(fig)


if __name__ == "__main__":
    fig_cases("fig1_manifold_cases.png")
    fig_scene_ts(
        "fig2_flat.png",
        "flat",
        [0.101, 0.103, 0.25, 0.9],
        (-0.8, 0.8),
        (-0.2, 0.55),
        "body1-seg(-3, 0, 3, 0)",
        tmax=0.6,
        title="flat landing (e = 0.5, μ = 0.3): both points have the same gap, they enter and leave the active set together",
    )
    fig_scene_ts(
        "fig3_tilted.png",
        "tilted",
        [0.05, 0.102, 0.256, 0.42, 1.2],
        (-0.8, 0.8),
        (-0.2, 0.75),
        "body1-seg(-3, 0, 3, 0)",
        tmax=0.8,
        title="corner-first landing, tilted 20° (e = 0.5, μ = 0.3): p0 (left corner) first, p1 (far corner) when the box has rocked down",
    )
    fig_scene_ts(
        "fig4_tipping.png",
        "tipping",
        [0.03, 0.047, 0.2, 0.3, 0.45],
        (-0.7, 1.2),
        (-0.8, 0.55),
        "body1-seg(-3, 0, 0, 0)",
        tmax=0.6,
        title="centre of mass past the segment end (e = 0, μ = 0.3): the box tips over the end vertex",
    )
    fig_static(
        "fig5_static.png",
        [
            ("overhang", 0.8, (-1.0, 0.4), (-0.25, 0.55), "overhang: CoM at x=-0.3"),
            (
                "junction",
                0.8,
                (-0.6, 0.8),
                (-0.25, 0.55),
                "two collinear segments, CoM at x=0.1",
            ),
            (
                "stack",
                0.8,
                (-1.2, 1.3),
                (-0.1, 0.8),
                "polytope on polytope: top box overhangs",
            ),
            (
                "corner",
                1.4,
                (-1.2, 1.2),
                (-0.1, 0.9),
                "square settled after corner landing",
            ),
        ],
    )
    fig_scene_ts(
        "fig6_corner.png",
        "corner",
        [0.08, 0.102, 0.25, 0.37, 1.0],
        (-0.6, 0.8),
        (0.3, 1.0),
        "body1-body2",
        tmax=1.0,
        title="polytope on polytope: a square (body 2) lands on its corner on a box (body 1), then topples flat",
    )
