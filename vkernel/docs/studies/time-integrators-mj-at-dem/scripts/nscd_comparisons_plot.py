# The tables and the plots of the comparisons of Moreau-Jean with the
# Anitescu-Tasora and explicit DEM schemes, from the CSV files written by
# nscd_comparisons (vkernel/soa/src/siconos/tests/nscd_comparisons.cpp).
# The study: vkernel/docs/studies/time-integrators-mj-at-dem.md; the cmake
# target vkernel_study_time_integrators runs both.
#
# usage: python3 nscd_comparisons_plot.py <directory of the CSV files>
# writes tables.md and, if matplotlib is there, one PNG per case, in the
# same directory. matplotlib is not in the build environment of siconos:
#   guix shell python python-matplotlib -- python3 nscd_comparisons_plot.py <dir>

import csv
import math
import os
import sys

directory = sys.argv[1]

# the schemes, in a fixed order: color (categorical slots 1-4 of the
# reference palette, validated for the light surface), marker
SCHEMES = ["Moreau-Jean", "AT (PGS)", "AT (APGD)", "DEM"]
COLORS = dict(zip(SCHEMES, ["#2a78d6", "#eb6834", "#1baf7a", "#eda100"]))
MARKERS = dict(zip(SCHEMES, ["o", "s", "^", "D"]))
SURFACE = "#fcfcfb"
INK = "#2b2b2b"
MUTED = "#8a8a86"


def read(name):
    path = os.path.join(directory, name + ".csv")
    if not os.path.exists(path):
        return []
    with open(path) as f:
        return list(csv.DictReader(f))


def num(x):
    return float(x)


def fmt(x, digits=2):
    x = float(x)
    if math.isinf(x) or math.isnan(x):
        return "unstable"
    return f"{x:.{digits}e}"


def table(title, header, rows):
    lines = [f"### {title}", "", "| " + " | ".join(header) + " |"]
    lines.append("|" + "|".join("---" for _ in header) + "|")
    lines += ["| " + " | ".join(str(c) for c in row) + " |" for row in rows]
    return "\n".join(lines) + "\n"


def by_scheme(rows, key):
    """rows grouped by scheme, then by key"""
    out = {}
    for row in rows:
        out.setdefault(row["scheme"], {})[row[key]] = row
    return out


md = ["# Moreau-Jean, Anitescu-Tasora, DEM: the comparisons", ""]

stability = read("stability")
if stability:
    g = by_scheme(stability, "h")
    hs = sorted({r["h"] for r in stability}, key=float)
    md.append(
        table(
            "Stability: a column of 5 disks at rest, 1 s - largest kinetic energy",
            ["h"] + SCHEMES,
            [
                [h]
                + [
                    fmt(g[s][h]["max_kinetic_energy"])
                    if g[s][h]["finite"] == "true"
                    else "unstable"
                    for s in SCHEMES
                ]
                for h in hs
            ],
        )
    )
    md.append(
        table(
            "Stability: largest overlap of a contact",
            ["h"] + SCHEMES,
            [
                [h]
                + [
                    fmt(g[s][h]["max_overlap"])
                    if g[s][h]["finite"] == "true"
                    else "unstable"
                    for s in SCHEMES
                ]
                for h in hs
            ],
        )
    )

columns = read("column")
for start in ("rest", "landing"):
    rows = [r for r in columns if r["start"] == start]
    if not rows:
        continue
    g = by_scheme(rows, "n")
    ns = sorted({r["n"] for r in rows}, key=int)
    what = "at rest" if start == "rest" else "landing from gaps of 0.01"
    md.append(
        table(
            f"Column of N disks {what}, 1 s - bottom overlap at the end (DEM: N m g / k_n)",
            ["N"] + SCHEMES + ["N m g / k_n"],
            [
                [n]
                + [fmt(g[s][n]["bottom_overlap_end"]) for s in SCHEMES]
                + [fmt(g["DEM"][n]["n_m_g_over_kn"])]
                for n in ns
            ],
        )
    )
    md.append(
        table(
            f"Column of N disks {what} - largest overlap of a contact, wall time (s)",
            ["N"] + [f"{s} overlap" for s in SCHEMES] + [f"{s} time" for s in SCHEMES],
            [
                [n]
                + [fmt(g[s][n]["max_overlap"]) for s in SCHEMES]
                + [f"{num(g[s][n]['wall_s']):.2f}" for s in SCHEMES]
                for n in ns
            ],
        )
    )

stick = read("incline_stick")
if stick:
    rows = []
    for s in SCHEMES:
        xs = [num(r["x"]) for r in stick if r["scheme"] == s]
        if xs:
            rows.append([s, fmt(max(abs(x) for x in xs)), fmt(xs[-1])])
    md.append(
        table(
            "Incline 20 deg, mu 0.5 (below the friction angle): displacement",
            ["scheme", "max abs(x)", "x(T)"],
            rows,
        )
    )


def error_table(name, title, value, exact):
    rows = read(name)
    if not rows:
        return
    g = by_scheme(rows, "h")
    hs = sorted({r["h"] for r in rows}, key=float, reverse=True)

    def err(r):
        if r["finite"] != "true":
            return "unstable"
        return fmt(num(r[value]) - num(r[exact]))

    md.append(
        table(
            f"{title} - error",
            ["h"] + SCHEMES,
            [[h] + [err(g[s][h]) for s in SCHEMES] for h in hs],
        )
    )
    md.append(
        table(
            f"{title} - largest lift-off",
            ["h"] + SCHEMES,
            [
                [h]
                + [
                    fmt(g[s][h]["max_lift"])
                    if g[s][h]["finite"] == "true"
                    else "unstable"
                    for s in SCHEMES
                ]
                for h in hs
            ],
        )
    )
    md.append(
        table(
            f"{title} - cost: wall time (s), solver iterations",
            ["h"] + SCHEMES,
            [
                [h]
                + [
                    f"{num(g[s][h]['wall_s']):.2f}, {g[s][h]['iterations']}"
                    for s in SCHEMES
                ]
                for h in hs
            ],
        )
    )


error_table(
    "sliding",
    "Sliding square, v0 = 2, mu = 0.3: stopping distance",
    "distance",
    "distance_exact",
)
error_table(
    "incline_slide",
    "Incline 30 deg, mu = 0.3 (above the friction angle): x(T = 1)",
    "x_end",
    "x_exact",
)

shear_summary = read("shear_summary")
if shear_summary:
    rows = []
    for mu in ("0", "0.3"):
        for s in SCHEMES:
            r = next(
                (r for r in shear_summary if r["scheme"] == s and r["mu"] == mu), None
            )
            if r:
                rows.append(
                    [
                        mu,
                        s,
                        r["h"],
                        fmt(r["max_above"]),
                        fmt(r["max_below"]),
                        r["iterations"],
                        f"{num(r['wall_s']):.2f}",
                    ]
                )
    md.append(
        table(
            "Dilatancy: a disk sheared over a fixed layer - height above / below the geometric one",
            [
                "mu",
                "scheme",
                "h",
                "max above",
                "max below",
                "solver iterations",
                "wall time (s)",
            ],
            rows,
        )
    )
    md.append(
        "Below the geometric path, Moreau-Jean and AT at mu = 0 (shear_detail.png):\n"
        "they keep u_N = 0, never the gap. In a hollow the path has a corner, the\n"
        "vertical velocity flips from -v tan 30 to +v tan 30 in one step; with\n"
        "theta = 1/2 the position moves with the mean of the two, about 0: the\n"
        "centre stays flat over that step while the path goes up, about 1e-3 below\n"
        "it (at most h v tan 30). Small teeth: the tangent motion over a convex arc\n"
        "opens the gap by ~(h v_T)^2 / (2 R) per step; when it is open, one step of\n"
        "free fall (g h^2), caught at the velocity level, the overlap kept.\n"
    )

dem_incline = read("dem_incline")
if dem_incline:
    md.append(
        table(
            "DEM on the incline 30 deg, mu 0.3: x(T = 1) ahead of the exact motion, vs k_t",
            ["rotation", "k_t", "excess", "loading estimate"],
            [
                [
                    r["rotation"],
                    fmt(r["kt"]),
                    fmt(r["excess"]),
                    fmt(r["loading_estimate"]),
                ]
                for r in dem_incline
            ],
        )
    )

with open(os.path.join(directory, "tables.md"), "w") as f:
    f.write("\n".join(md))
print(f"wrote {os.path.join(directory, 'tables.md')}")

try:
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.patches
    import matplotlib.pyplot as plt
except ImportError:
    print("no matplotlib: tables only")
    sys.exit(0)

plt.rcParams.update(
    {
        "figure.facecolor": SURFACE,
        "axes.facecolor": SURFACE,
        "axes.edgecolor": MUTED,
        "axes.labelcolor": INK,
        "text.color": INK,
        "xtick.color": MUTED,
        "ytick.color": MUTED,
        "axes.grid": True,
        "grid.color": "#e6e6e3",
        "grid.linewidth": 0.6,
        "axes.spines.top": False,
        "axes.spines.right": False,
        "legend.frameon": False,
        "lines.linewidth": 1.5,
        "lines.markersize": 6,
    }
)


def schemes_chart(ax, data, floor=None, every=None):
    """data[scheme] = (x, y): one line per scheme, Moreau-Jean on top; the
    markers staggered so that lines with the same values (to 1e-6) stay
    visible, and named in a footnote under the chart; on a log axis, the
    values below floor drawn at floor, and named"""
    notes = []
    names = [s for s in SCHEMES if s in data]
    for k, s in enumerate(names):
        x, y = data[s]
        if floor is not None:
            low = [v for v in y if v < floor]
            if low and len(low) == len(y) and all(v <= 0 for v in low):
                notes.append(f"{s}: 0, drawn at {floor:.0e}")
            elif low:
                zero = " (0)" if all(v <= 0 for v in low) else ""
                notes.append(
                    f"{s}: the values below {floor:.0e}{zero} drawn at {floor:.0e}"
                )
            y = [max(v, floor) for v in y]
        n = len(x)
        step = every or max(1, n)  # without every: a marker on every point
        mark = (k * max(1, step // len(names)) % step, step) if every else None
        ax.plot(
            x,
            y,
            color=COLORS[s],
            marker=MARKERS[s],
            label=s,
            markeredgecolor=SURFACE,
            markeredgewidth=0.8,
            markevery=mark,
            zorder=10 - k,
        )
    # the series with the same values (to 1e-9 relative)
    same = []
    for i, a in enumerate(names):
        for b in names[i + 1 :]:
            ya, yb = data[a][1], data[b][1]
            if len(ya) == len(yb) and all(
                abs(u - v) <= 1e-6 * max(abs(u), abs(v), 1e-300) for u, v in zip(ya, yb)
            ):
                same.append((a, b))
    groups = []
    for a, b in same:
        for g in groups:
            if a in g or b in g:
                g.update((a, b))
                break
        else:
            groups.append({a, b})
    for g in groups:
        notes.append(
            "same values (to 1e-6): " + ", ".join(s for s in SCHEMES if s in g)
        )
    if notes:
        # under the x label, outside the plot (saved with bbox_inches="tight")
        ax.annotate(
            "\n".join(notes),
            xy=(0, 0),
            xycoords="axes fraction",
            xytext=(0, -38),
            textcoords="offset points",
            fontsize=7,
            color=MUTED,
            va="top",
            ha="left",
            annotation_clip=False,
        )


def save(fig, name):
    fig.tight_layout()
    path = os.path.join(directory, name + ".png")
    fig.savefig(path, dpi=150, bbox_inches="tight", facecolor=SURFACE)
    plt.close(fig)
    print(f"wrote {path}")


# the scenes: a schematic of each case, in its own panel (bodies, ground,
# gravity, the imposed motion, what the schemes do)

BODY = "#e9e9e6"


def scene_axes(ax, title):
    ax.set_aspect("equal")
    ax.axis("off")
    ax.set_title(title, loc="left")


def rotate(x, y, angle):
    c, s_ = math.cos(angle), math.sin(angle)
    return c * x - s_ * y, s_ * x + c * y


def ground(ax, x0, x1, angle=0.0):
    """the ground from x0 to x1 along a line through the origin at angle"""
    (ax0, ay0), (ax1, ay1) = rotate(x0, 0, angle), rotate(x1, 0, angle)
    ax.plot([ax0, ax1], [ay0, ay1], color=INK, linewidth=1.5, solid_capstyle="butt")
    n = int((x1 - x0) / 0.35)
    for k in range(n + 1):
        u = x0 + k * (x1 - x0) / max(n, 1)
        bx, by = rotate(u, 0, angle)
        tx, ty = rotate(u - 0.2, -0.25, angle)
        ax.plot([bx, tx], [by, ty], color=MUTED, linewidth=0.8)


def square(ax, u, v, angle=0.0, ls="-", alpha=1.0, edge=INK):
    """a unit square resting with its base at height v, centred at u along
    the ground at angle"""
    corners = [
        rotate(u + dx, v + dy, angle)
        for dx, dy in ((-0.5, 0), (0.5, 0), (0.5, 1), (-0.5, 1))
    ]
    ax.add_patch(
        matplotlib.patches.Polygon(
            corners,
            closed=True,
            facecolor=BODY if ls == "-" else "none",
            edgecolor=edge,
            linestyle=ls,
            linewidth=1.2,
            alpha=alpha,
        )
    )
    return rotate(u, v + 0.5, angle)


def disk(ax, x, y, ls="-"):
    ax.add_patch(
        matplotlib.patches.Circle(
            (x, y),
            0.5,
            facecolor=BODY if ls == "-" else "none",
            edgecolor=INK,
            linestyle=ls,
            linewidth=1.2,
        )
    )


def arrow(ax, x, y, dx, dy, text="", color=INK, ls="-", where="end"):
    ax.annotate(
        "",
        xy=(x + dx, y + dy),
        xytext=(x, y),
        arrowprops=dict(
            arrowstyle="-|>",
            color=color,
            linewidth=1.3,
            linestyle=ls,
            shrinkA=0,
            shrinkB=0,
        ),
    )
    if text:
        tx, ty = (x + dx, y + dy) if where == "end" else (x + dx / 2, y + dy / 2)
        ax.text(tx + 0.08, ty, text, color=INK, fontsize=8, va="center")


def gravity(ax, x, y, length=1.0):
    arrow(ax, x, y, 0, -length, "g", color=MUTED)


def note(ax, x, y, text, **kw):
    ax.text(x, y, text, color=INK, fontsize=8, va="top", **kw)


def scene_stability(ax):
    scene_axes(ax, "The scene")
    ground(ax, -2, 2)
    for i in range(5):
        disk(ax, 0, 0.5 + i)
    gravity(ax, 1.4, 4.5)
    note(
        ax,
        -2.2,
        -0.5,
        "5 disks at rest, in contact; 1 s\nh from 1e-4 to 1e-2\nDEM: a penalty spring at each\ncontact, stable for h < ~2 sqrt(m/k_n)",
    )
    ax.set_xlim(-2.4, 2.6)
    ax.set_ylim(-2.4, 5.6)


def scene_column(ax):
    scene_axes(ax, "The scenes")
    ground(ax, -3.2, 3.2)
    for i in range(4):
        disk(ax, -1.6, 0.5 + i)
    gap = 0.3  # 0.01, drawn x30
    for i in range(4):
        disk(ax, 1.6, 0.5 + gap + (1 + gap) * i)
    ax.text(-1.6, 4.4, "at rest", ha="center", fontsize=8, color=INK)
    ax.text(1.6, 5.6, "landing", ha="center", fontsize=8, color=INK)
    gravity(ax, 3.0, 5.0)
    note(
        ax,
        -3.3,
        -0.5,
        "N disks (1 to 50), 1 s; gaps of 0.01\n(drawn x30). Moreau-Jean, AT: no\noverlap at rest; DEM: the bottom\ncontact carries N m g, overlap N m g / k_n",
    )
    ax.set_xlim(-3.4, 3.6)
    ax.set_ylim(-2.6, 6.2)


def scene_incline(ax, alpha_deg, mu, slides):
    a = math.radians(alpha_deg)
    scene_axes(ax, "The scene")
    ground(ax, -2.5, 2.5, -a)
    cx, cy = square(ax, -0.6, 0, -a)
    gravity(ax, cx, cy, 0.9)
    if slides:
        dx, dy = rotate(1.4, 0, -a)
        bx, by = rotate(-0.6, 1.25, -a)
        arrow(
            ax,
            bx,
            by,
            dx,
            dy,
            "a = g (sin \u03b1 \u2212 \u03bc cos \u03b1)",
            where="end",
        )
        square(ax, 1.4, 0, -a, ls="--", edge=MUTED)
        text = (
            f"\u03b1 = {alpha_deg}\u00b0: tan \u03b1 = {math.tan(a):.2f} > \u03bc = {mu}: slides\n"
            "Moreau-Jean: x = a t\u00b2 / 2 exactly\nAT: lifts off (u_N = \u03bc |u_T|)"
        )
    else:
        bx, by = rotate(-0.6, 1.25, -a)
        dx, dy = rotate(0.9, 0, -a)
        arrow(ax, bx, by, dx, dy, "x ?", color=MUTED, ls="--", where="end")
        text = (
            f"\u03b1 = {alpha_deg}\u00b0: tan \u03b1 = {math.tan(a):.2f} < \u03bc = {mu}: sticks\n"
            "Moreau-Jean, AT: x = 0\nDEM: the tangential springs stretch"
        )
    note(ax, -2.6, -1.6, text)
    ax.set_xlim(-2.8, 2.8)
    ax.set_ylim(-3.0, 2.4)


def scene_sliding(ax):
    scene_axes(ax, "The scene")
    ground(ax, -1.2, 3.6)
    cx, cy = square(ax, 0, 0)
    arrow(ax, 0.55, 0.5, 0.9, 0, "v\u2080 = 2")
    arrow(ax, 0.2, 0.12, -0.7, 0)
    ax.text(-0.95, 0.3, "\u03bc N", fontsize=8, color=INK)
    square(ax, 2.8, 0, ls="--", edge=MUTED)
    ax.annotate(
        "",
        xy=(2.8, -0.7),
        xytext=(0, -0.7),
        arrowprops=dict(arrowstyle="<->", color=MUTED, linewidth=0.9),
    )
    ax.text(
        1.4,
        -0.95,
        "d = v\u2080\u00b2 / (2 \u03bc g) = 0.68",
        ha="center",
        fontsize=8,
        color=INK,
    )
    # the hop of AT, schematic
    xs = [0.4 + 0.02 * k for k in range(60)]
    ax.plot(
        xs,
        [1.0 + 0.9 * math.sin(math.pi * (x - 0.4) / 1.2) for x in xs],
        color=COLORS["AT (PGS)"],
        linestyle=":",
        linewidth=1.3,
    )
    ax.text(1.0, 2.05, "AT: hops", ha="center", fontsize=8, color=INK)
    gravity(ax, -0.9, 2.2)
    note(
        ax,
        -1.3,
        -1.25,
        "a unit square launched at v\u2080 on the ground,\n\u03bc = 0.3: Moreau-Jean stays on the ground\nand stops at d",
    )
    ax.set_xlim(-1.4, 3.8)
    ax.set_ylim(-2.6, 2.6)


if stability:
    fig, (sc, ax) = plt.subplots(
        1, 2, figsize=(10, 4), gridspec_kw={"width_ratios": [1, 1.6]}
    )
    scene_stability(sc)
    g = by_scheme(stability, "h")
    data = {}
    for s in SCHEMES:
        hs = sorted(g[s], key=float)
        ok = [h for h in hs if g[s][h]["finite"] == "true"]
        data[s] = (
            [num(h) for h in ok],
            [num(g[s][h]["max_kinetic_energy"]) for h in ok],
        )
    schemes_chart(ax, data, floor=1e-30)
    for s in SCHEMES:
        hs = sorted(g[s], key=float)
        bad = [h for h in hs if g[s][h]["finite"] != "true"]
        if bad:
            ax.axvspan(
                num(bad[0]), num(hs[-1]) * 1.3, color=COLORS[s], alpha=0.08, linewidth=0
            )
            ax.text(
                num(bad[0]) * 1.05,
                1e-2,
                f"{s}: unstable",
                color=INK,
                fontsize=8,
                rotation=90,
                va="bottom",
            )
    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_xlabel("time step h")
    ax.set_ylabel("largest kinetic energy over 1 s (J)")
    ax.set_title("The column should stay at rest", loc="left")
    ax.legend()
    save(fig, "stability")

if columns:
    fig, axes = plt.subplots(
        1, 3, figsize=(14, 4), gridspec_kw={"width_ratios": [1, 1.3, 1.3]}
    )
    scene_column(axes[0])
    for ax, start in zip(axes[1:], ("rest", "landing")):
        rows = [r for r in columns if r["start"] == start]
        g = by_scheme(rows, "n")
        data = {}
        for s in SCHEMES:
            ns = sorted(g[s], key=int)
            data[s] = ([int(n) for n in ns], [num(g[s][n]["max_overlap"]) for n in ns])
        schemes_chart(ax, data, floor=1e-16)
        ns = sorted(g["DEM"], key=int)
        ax.plot(
            [int(n) for n in ns],
            [num(g["DEM"][n]["n_m_g_over_kn"]) for n in ns],
            color=MUTED,
            linestyle="--",
            linewidth=1,
            label="N m g / k_n",
        )
        ax.set_xscale("log")
        ax.set_yscale("log")
        ax.set_xlabel("number of disks N")
        ax.set_ylabel("largest overlap of a contact (m)")
        ax.set_title(
            "at rest" if start == "rest" else "landing from gaps of 0.01", loc="left"
        )
    axes[1].legend()
    save(fig, "column")

if stick:
    fig, (sc, ax) = plt.subplots(
        1, 2, figsize=(10, 4), gridspec_kw={"width_ratios": [1, 1.6]}
    )
    scene_incline(sc, 20, 0.5, False)
    data = {}
    for s in SCHEMES:
        rows = [r for r in stick if r["scheme"] == s]
        data[s] = ([num(r["t"]) for r in rows], [abs(num(r["x"])) for r in rows])
    schemes_chart(ax, data, floor=1e-16, every=20)
    ax.set_yscale("log")
    ax.set_xlabel("time (s)")
    ax.set_ylabel("abs(displacement down the incline) (m)")
    ax.set_title("The square should not move", loc="left")
    ax.legend()
    save(fig, "incline_stick")

for name, title, value, exact, draw in (
    (
        "sliding",
        "Sliding square: stopping distance",
        "distance",
        "distance_exact",
        scene_sliding,
    ),
    (
        "incline_slide",
        "Incline 30 deg, mu 0.3: x(T = 1)",
        "x_end",
        "x_exact",
        lambda ax: scene_incline(ax, 30, 0.3, True),
    ),
):
    rows = read(name)
    paths = read(name + "_path")
    if paths:
        # the scene and the trajectories: the lift-off along the way
        fig, (sc, ax) = plt.subplots(
            1, 2, figsize=(10, 4), gridspec_kw={"width_ratios": [1, 1.6]}
        )
        draw(sc)
        data = {}
        for s in SCHEMES:
            p = [r for r in paths if r["scheme"] == s]
            data[s] = ([num(r["x"]) for r in p], [num(r["lift"]) for r in p])
        schemes_chart(ax, data, every=24)
        ax.set_xlabel(
            "distance travelled x (m)"
            if name == "sliding"
            else "distance down the incline x (m)"
        )
        ax.set_ylabel("lift-off above the rest height (m)")
        ax.set_title("The path of the square (h = 5e-3, DEM 1e-4)", loc="left")
        ax.legend()
        save(fig, name)
    if not rows:
        continue
    g = by_scheme(rows, "h")
    fig, axes = plt.subplots(1, 2, figsize=(11, 4))
    for ax, ylabel, get, what in (
        (
            axes[0],
            "abs(error) (m)",
            lambda r: abs(num(r[value]) - num(r[exact])),
            "error",
        ),
        (axes[1], "largest lift-off (m)", lambda r: num(r["max_lift"]), "lift-off"),
    ):
        data = {}
        for s in SCHEMES:
            hs = sorted((h for h in g[s] if g[s][h]["finite"] == "true"), key=float)
            data[s] = ([num(h) for h in hs], [get(g[s][h]) for h in hs])
        schemes_chart(ax, data, floor=1e-16)
        ax.set_xscale("log")
        ax.set_yscale("log")
        ax.set_xlabel("time step h")
        ax.set_ylabel(ylabel)
        ax.set_title(f"{title}: {what}", loc="left")
    axes[0].legend()
    save(fig, name + "_h")


def scene_shear(ax):
    scene_axes(ax, "The scene")
    ground(ax, -0.6, 5.6)
    for i in range(6):
        disk(ax, i * 1.0, 0.5)
        ax.add_patch(
            matplotlib.patches.Circle(
                (i * 1.0, 0.55),
                0.5,
                facecolor="none",
                edgecolor=MUTED,
                hatch="////",
                linewidth=0,
            )
        )
    x0 = 1.5
    disk(ax, x0, 0.5 + math.sqrt(0.75))
    arrow(ax, x0 + 0.55, 0.5 + math.sqrt(0.75), 0.8, 0, "v imposed")
    # the geometric path of its centre
    xs = [x0 + 0.02 * k for k in range(151)]
    ys = [0.5 + max(math.sqrt(max(1 - (x - i) ** 2, 0)) for i in range(6)) for x in xs]
    ax.plot(xs, ys, color=MUTED, linestyle="--", linewidth=1)
    ax.text(
        3.0, 2.35, "geometric path of the centre", fontsize=8, color=INK, ha="center"
    )
    gravity(ax, -0.4, 3.4)
    note(
        ax,
        -0.7,
        -0.5,
        "a disk sheared over a layer of fixed disks\n(hatched), rotation locked: its contacts slide;\n"
        "from a hollow (\u221a3 r) to a crest (2 r): the\ngeometric dilatancy. AT, \u03bc > 0: u_N = \u03bc |u_T|\n"
        "launches it at each sliding step, it falls\nback at the next: above the path",
    )
    ax.set_xlim(-0.8, 5.8)
    ax.set_ylim(-2.4, 3.6)


shear_rows = read("shear")
if shear_rows:
    fig, axes = plt.subplots(
        1, 3, figsize=(15, 4), gridspec_kw={"width_ratios": [1.1, 1.3, 1.3]}
    )
    scene_shear(axes[0])
    for ax, mu in zip(axes[1:], ("0", "0.3")):
        data = {}
        for s in SCHEMES:
            rows = [r for r in shear_rows if r["scheme"] == s and r["mu"] == mu]
            data[s] = (
                [num(r["x"]) for r in rows],
                [num(r["y"]) - num(r["y_geometric"]) for r in rows],
            )
        schemes_chart(ax, data, every=40)
        ax.axhline(0, color=MUTED, linewidth=0.8, linestyle="--")
        ax.set_xlabel("position of the sheared disk x (m)")
        ax.set_ylabel("height above the geometric one (m)")
        ax.set_title(f"\u03bc = {mu}", loc="left")
    axes[1].legend()
    save(fig, "shear")


def shear_geometric(x):
    """the height of the centre of a disk resting on the layer (r = 0.5,
    fixed disks at the integers, at the height r), as in nscd_comparisons"""
    return 0.5 + max(math.sqrt(1 - (x - i) ** 2) for i in range(20) if abs(x - i) < 1)


# the details of Moreau-Jean, mu = 0: the trajectory itself, a hollow and a
# small tooth, one point per step
mj0 = [r for r in shear_rows if r["scheme"] == "Moreau-Jean" and r["mu"] == "0"]
if mj0:
    xs = [num(r["x"]) for r in mj0]
    ys = [num(r["y"]) for r in mj0]
    fig, axes = plt.subplots(1, 3, figsize=(15, 4))
    color = COLORS["Moreau-Jean"]

    def geometric_line(ax, x0, x1):
        gx = [x0 + (x1 - x0) * k / 400 for k in range(401)]
        # the corners of the hollows, exactly
        gx = sorted(set(gx) | {i + 0.5 for i in range(20) if x0 < i + 0.5 < x1})
        ax.plot(
            gx,
            [shear_geometric(x) for x in gx],
            color=MUTED,
            linestyle="--",
            linewidth=1,
            label="geometric path",
        )

    ax = axes[0]
    geometric_line(ax, xs[0], xs[-1])
    ax.plot(xs, ys, color=color)
    ax.set_ylabel("height of the centre y (m)")
    ax.set_title("\u03bc = 0: Moreau-Jean, the centre (dashed: geometric)", loc="left")

    ax = axes[1]
    sel = [i for i, x in enumerate(xs) if 6.48 < x < 6.53]
    geometric_line(ax, 6.48, 6.53)
    ax.plot(
        [xs[i] for i in sel],
        [ys[i] for i in sel],
        color=color,
        marker="o",
        markersize=4,
        markeredgecolor=SURFACE,
        label="centre, one point per step",
    )
    ax.set_ylabel("height of the centre y (m)")
    ax.set_title("the hollow x = 6.5: the step over the corner", loc="left")
    ax.legend(loc="upper left")
    note(
        ax,
        6.5085,
        1.3692,
        "v_y: -v tan 30 \u2192 +v tan 30\nin one step; \u03b8 = 1/2: moves with\n"
        "the mean, \u2248 0: flat, then ~1e-3\nbelow the path, kept:\n"
        "u_N = 0, not the gap",
    )

    ax = axes[2]
    sel = [i for i, x in enumerate(xs) if 5.68 < x < 5.83]
    ax.plot(
        [xs[i] for i in sel],
        [ys[i] - shear_geometric(xs[i]) for i in sel],
        color=color,
        marker="o",
        markersize=2,
        markeredgewidth=0,
    )
    ax.axhline(0, color=MUTED, linewidth=0.8, linestyle="--")
    ax.set_ylabel("height above the geometric one (m)")
    ax.set_title("small teeth: the gap opens, one step of free fall", loc="left")
    note(
        ax,
        5.685,
        -2.4e-4,
        "tangent motion on a convex arc: the gap\nopens, ~(h v_T)\u00b2/(2R) per step; open:\n"
        "one step of free fall (g h\u00b2), caught\nat the velocity level",
    )
    ax.set_ylim(-3.4e-4, 3e-5)
    for ax in axes:
        ax.set_xlabel("position of the sheared disk x (m)")
    save(fig, "shear_detail")

if dem_incline:
    fig, ax = plt.subplots(figsize=(6.4, 4))
    for rotation, marker, ls in (("free", "D", "-"), ("locked", "o", ":")):
        rows = [r for r in dem_incline if r["rotation"] == rotation]
        ax.plot(
            [num(r["kt"]) for r in rows],
            [num(r["excess"]) for r in rows],
            color=COLORS["DEM"],
            marker=marker,
            linestyle=ls,
            markeredgecolor=SURFACE,
            label=f"DEM, rotation {rotation}",
        )
    rows = [r for r in dem_incline if r["rotation"] == "free"]
    ax.plot(
        [num(r["kt"]) for r in rows],
        [num(r["loading_estimate"]) for r in rows],
        color=MUTED,
        linestyle="--",
        linewidth=1,
        label="loading of the tangential springs (estimate)",
    )
    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_xlabel("tangential stiffness k_t (N/m)")
    ax.set_ylabel("x(T = 1) ahead of the exact motion (m)")
    ax.set_title("DEM on the incline: the start of the springs", loc="left")
    ax.legend()
    save(fig, "dem_incline")
