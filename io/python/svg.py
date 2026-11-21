# Siconos is a program dedicated to modeling, simulation and control
# of non smooth dynamical systems.
#
# Copyright 2026 INRIA.
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
# http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.
#
"""SVG pictures of a 2D mechanics simulation.

One picture is one output step: the static objects, the bodies with their
shape coloured by their speed, the contacts (the reaction impulse on the
first body, the Coulomb cone about the normal). The fixed segments are
hatched on their material side, as a floor in a drawing.

The geometry is the one of the ParaView export (vtkhdf.ShapeTable): the
shapes, the radii, the times of birth and death of the static objects.

A segment has no side in the model. Its material side is inferred from
the contacts with it: the normal of a contact points from the obstacle
into the body, so the material is on the other side. A contact with a
fixed obstacle has ds1 == ds2 (the body) and point_a on the obstacle: it
is matched to the segment by distance. The votes are those of the output
steps up to the drawn one; a segment without contact takes the side
opposite to the bodies near it, at the drawn step then at the first one;
a segment with votes on both sides, or none, is not hatched and a warning
names it. material_sides overrides the inference.

Two callers write the same bytes for the same step:
- convert(io_filename, ...) from the HDF5 file, siconos_svg on the
  command line;
- the svg_output run option of MechanicsHdf5Runner, during the run, from
  the arrays of the IO, rounded to float32 as the HDF5 tables.

The file depends on numpy (and h5py to read a file). A PNG is a
conversion: rsvg-convert -o out.png out.svg (guix shell librsvg).
"""

import os
import sys

import numpy as np
from siconos.io import vtkhdf

# the cf columns, after the time column
CF = vtkhdf._CF

INK = "#2b2b2b"
MUTED = "#8c8c8c"
BACKGROUND = "#fcfcfb"
STATIC_FILL = "#d9d9d6"
FORCE = "#c0392b"
CONE = "#e8a33d"
# a fixed ramp (viridis stops) for the speed of the bodies; the legend gives
# the range of the drawn step
RAMP = ("#440154", "#3b528b", "#21918c", "#5ec962", "#fde725")

WIDTH = 800.0
TITLE_HEIGHT = 22

# the body modes: one element per body below CIRCLE_LIMIT bodies (a circle
# and a rotation tick for a disk), one <use> of a shape in <defs> per body
# below USE_LIMIT, above the bodies merged in one path per colour
CIRCLE_LIMIT = 2000
USE_LIMIT = 20000
# the contact arrows are drawn below ARROW_LIMIT contacts, the cones below
# CONE_LIMIT, unless asked
ARROW_LIMIT = 2000
CONE_LIMIT = 20

# the hatches: a stroke every HATCH_STEP typical body sizes, at least one
# every 1/HATCH_MAX of the view diagonal
HATCH_STEP = 0.5
HATCH_MAX = 150
# a contact point closer to a segment than this fraction of the scene
# diagonal is on it (the tables are float32)
ON_SEGMENT = 1e-5
# a contact normal closer to the segment direction than this cosine is at an
# end point: no vote
END_POINT_COSINE = 0.1


def scheme_name(osi):
    """The name of an integrator (a class or an object) for the title, as
    the run options of the file have it."""
    if osi is None or isinstance(osi, str):
        return osi
    return getattr(osi, "__name__", type(osi).__name__)


def _file_scheme(data):
    """The integrator of the run options of a file, None if not written."""
    import json

    options = data.get("siconos_mechanics_run_options")
    if options is None or "options" not in options.attrs:
        return None
    return json.loads(options.attrs["options"]).get("osi")


def warn(message):
    print(f"siconos svg: warning: {message}", file=sys.stderr)


def svg_filename(io_filename, time, directory=None):
    """<stem>_svg_<time>.svg, beside the HDF5 file or in directory. The time
    is written with 6 decimals, zero padded so that the names sort."""
    stem = os.path.splitext(os.path.basename(io_filename))[0]
    directory = directory if directory is not None else os.path.dirname(io_filename)
    return os.path.join(directory, f"{stem}_svg_{time:013.6f}.svg")


def parse_times(text):
    """The requests of svg_output: "0.5,1,2.0" -> [0.5, 1, 2.0]. A value with
    a decimal point or an exponent is a time (float), else a simulation
    step number (int)."""
    if isinstance(text, (list, tuple)):
        return list(text)
    requests = []
    for item in str(text).strip("[] ").split(","):
        item = item.strip()
        if not item:
            continue
        if any(c in item for c in ".eE"):
            requests.append(float(item))
        else:
            requests.append(int(item))
    return requests


def as_table(values, columns):
    """values rounded as in the HDF5 tables (float32), as float64."""
    values = np.asarray(values, dtype=np.float32).astype(np.float64)
    return values.reshape(-1, columns)


# -- the material side of the fixed segments ---------------------------------


def _segments(geometry, ids):
    """(key, p1, p2) of the two-point cells of a static geometry, key =
    (static id, rank of the cell in its object)."""
    segments, rank = [], {}
    for cell, i in zip(geometry.cells, ids):
        i = int(i)
        k = rank.get(i, 0)
        rank[i] = k + 1
        if len(cell) == 2:
            p = geometry.points[np.asarray(cell)][:, :2]
            segments.append(((i, k), p[0], p[1]))
    return segments


def _distance_to_segment(points, p1, p2):
    """The distances of points to [p1, p2] and the parameters of their
    projections."""
    d = p2 - p1
    t = np.clip(((points - p1) @ d) / (d @ d), 0.0, 1.0)
    return np.linalg.norm(points - (p1 + t[:, None] * d), axis=1), t


def _perp(p1, p2):
    d = p2 - p1
    return np.array((-d[1], d[0])) / np.linalg.norm(d)


class Sides:
    """The votes for the side of the bodies of each fixed segment: +1 on the
    side of perp(p2 - p1), -1 on the other."""

    def __init__(self, tolerance, forced=None):
        self.tolerance = tolerance
        # key -> [votes +1, votes -1]
        self.votes = {}
        # static id -> +1 or -1: the material side, given
        self.forced = dict(forced or {})

    def add_contacts(self, segments, contacts):
        """The votes of the contacts (cf rows without time) of a step."""
        fixed = contacts[contacts[:, CF["ds1"]] == contacts[:, CF["ds2"]]]
        if len(fixed) == 0:
            return
        points = fixed[:, CF["point_a"]][:, :2]
        normals = fixed[:, CF["normal"]][:, :2]
        for key, p1, p2 in segments:
            distance, _ = _distance_to_segment(points, p1, p2)
            dots = normals[distance < self.tolerance] @ _perp(p1, p2)
            dots = dots[np.abs(dots) > END_POINT_COSINE]
            if len(dots):
                count = self.votes.setdefault(key, [0, 0])
                count[0] += int(np.sum(dots > 0))
                count[1] += int(np.sum(dots < 0))

    def material(self, segments, centre_sets, size):
        """key -> unit normal pointing into the material, or None (not
        hatched), with a warning for each segment not hatched.

        centre_sets: the body centres to try, in order, for a segment
        without contact.
        """
        sides = {}
        for key, p1, p2 in segments:
            perp = _perp(p1, p2)
            if key[0] in self.forced:
                sides[key] = self.forced[key[0]] * perp
                continue
            plus, minus = self.votes.get(key, (0, 0))
            source = "contacts"
            if plus + minus == 0:
                source = "bodies near it"
                for centres in centre_sets:
                    distance, t = _distance_to_segment(centres, p1, p2)
                    near = (distance < 2 * size) & (t > 0) & (t < 1)
                    side = (centres[near] - p1) @ perp
                    plus, minus = int(np.sum(side > 0)), int(np.sum(side < 0))
                    if plus + minus:
                        break
            where = f"static {key[0]} segment {_point(p1)} {_point(p2)}"
            if plus + minus == 0:
                warn(f"{where}: no contact nor body near it, not hatched")
                sides[key] = None
            elif plus and minus:
                warn(
                    f"{where}: bodies on both sides ({source}: {plus} and "
                    f"{minus}), not hatched"
                )
                sides[key] = None
            else:
                # the material is on the other side of the bodies
                sides[key] = -perp if plus else perp
        return sides


def _point(p):
    return f"({p[0]:.6g}, {p[1]:.6g})"


# -- the picture ----------------------------------------------------------------


class _View:
    """World (x, y) to SVG px, y up, 2 decimals."""

    def __init__(self, low, high, margin):
        self.x0, self.y1 = low[0] - margin, high[1] + margin
        width = high[0] - low[0] + 2 * margin
        height = high[1] - low[1] + 2 * margin
        self.scale = WIDTH / width
        self.width, self.height = WIDTH, height * self.scale

    def x(self, p):
        return f"{(p[0] - self.x0) * self.scale:.2f}"

    def y(self, p):
        return f"{(self.y1 - p[1]) * self.scale:.2f}"

    def xy(self, p):
        return f"{self.x(p)},{self.y(p)}"

    def points(self, ps):
        return " ".join(self.xy(p) for p in ps)

    def length(self, v):
        return f"{v * self.scale:.2f}"

    def delta(self, d):
        return f"{d[0] * self.scale:.2f},{-d[1] * self.scale:.2f}"


def _rotate(angles, points):
    """points (N, 2) rotated by angles (N,) or a scalar."""
    c, s = np.cos(angles), np.sin(angles)
    x, y = points[..., 0], points[..., 1]
    return np.stack((c * x - s * y, s * x + c * y), axis=-1)


def _ramp_index(values, low, high):
    if high <= low:
        return np.zeros(len(values), dtype=np.int64)
    u = (values - low) / (high - low) * (len(RAMP) - 1)
    return np.clip(np.rint(u), 0, len(RAMP) - 1).astype(np.int64)


def _escape(text):
    return text.replace("&", "&amp;").replace("<", "&lt;").replace(">", "&gt;")


class Picture:
    """The pictures of a simulation: its shapes (a vtkhdf.ShapeTable), the
    votes for the sides of its segments, the drawing options.

    mode: "auto", "circle", "use" or "bins" (see CIRCLE_LIMIT).
    arrows, cones: True, False or None (by the number of contacts).
    material_sides: {static id: +1 | -1} the material side of the
    segments of a static object, on the side of perp(p2 - p1) (+1) or
    the other, instead of the inference.
    scheme: the name of the integrator, in the title.
    """

    def __init__(
        self,
        shapes,
        mode="auto",
        arrows=None,
        cones=None,
        material_sides=None,
        labels=False,
        scheme=None,
    ):
        if shapes.dimension != 2:
            raise ValueError("SVG pictures are for 2D simulations")
        self.shapes = shapes
        self.mode = mode
        self.arrows = arrows
        self.cones = cones
        self.labels = labels
        self.scheme = scheme
        bounds = shapes.bounds or [0, 1, 0, 1, 0, 0]
        self._diagonal = float(np.hypot(bounds[1] - bounds[0], bounds[3] - bounds[2]))
        self.sides = Sides(ON_SEGMENT * self._diagonal, material_sides)
        self._unit = shapes._block_of_key.get(("unit",))
        self._first_centres = None

    def segments(self, static_rows, time):
        geometry, ids = self.shapes.static_geometry(static_rows, time)
        return geometry, ids, _segments(geometry, ids)

    def vote(self, time, static_rows, contacts):
        """Count the contacts of an output step in the votes."""
        _, _, segments = self.segments(static_rows, time)
        if len(contacts):
            self.sides.add_contacts(segments, contacts)

    def first_bodies(self, bodies):
        """The bodies of the first output step, for the segments without
        contact."""
        if self._first_centres is None:
            self._first_centres = np.asarray(bodies)[:, 1:3].copy()

    def write(
        self, out, time, bodies, velocities, static_rows, contacts, requested=None
    ):
        """Write the picture of an output step in the text file out.

        bodies: (N, 8) rows id, x, y, z, q0, q1, q2, q3 (the dynamic table
        without its time column); velocities: (N, 7) id, v, w; static_rows:
        (M, 8) as bodies; contacts: the cf rows without time. The votes of
        the step must have been counted (vote).
        """
        self._largest = 0.0
        order = np.argsort(bodies[:, 0], kind="stable")
        bodies = bodies[order]
        velocities = velocities[np.argsort(velocities[:, 0], kind="stable")]
        ids = bodies[:, 0].astype(np.int64)
        centres = bodies[:, 1:3]
        angles = 2 * np.arctan2(bodies[:, 7], bodies[:, 4])
        data = self.shapes.point_data(ids)
        size = self.shapes.typical_size

        geometry, cell_ids, segments = self.segments(static_rows, time)
        centre_sets = [centres] + (
            [self._first_centres] if self._first_centres is not None else []
        )
        sides = self.sides.material(segments, centre_sets, size)

        view = self._view(centres, data["disk_radius"], geometry, size)
        mode = self.mode
        if mode == "auto":
            n = len(ids)
            mode = (
                "circle" if n <= CIRCLE_LIMIT else "use" if n <= USE_LIMIT else "bins"
            )
        arrows = len(contacts) <= ARROW_LIMIT if self.arrows is None else self.arrows
        cones = len(contacts) <= CONE_LIMIT if self.cones is None else self.cones

        speed = np.linalg.norm(velocities[:, 1:3], axis=1)
        low, high = (
            (float(speed.min()), float(speed.max())) if len(speed) else (0.0, 0.0)
        )
        colours = _ramp_index(speed, low, high)

        w = out.write
        height = view.height + TITLE_HEIGHT
        w('<?xml version="1.0" encoding="UTF-8"?>\n')
        w(
            f'<svg xmlns="http://www.w3.org/2000/svg" '
            f'xmlns:xlink="http://www.w3.org/1999/xlink" '
            f'width="{view.width:.0f}" height="{height:.0f}" '
            f'viewBox="0 {-TITLE_HEIGHT} {view.width:.0f} {height:.0f}" '
            f'font-family="sans-serif">\n'
        )
        w(
            f'<rect x="0" y="{-TITLE_HEIGHT}" width="100%" height="100%" fill="{BACKGROUND}"/>\n'
        )
        self._write_static(w, view, geometry, cell_ids, sides, size)
        self._write_bodies(w, view, mode, ids, centres, angles, data, colours)
        if len(contacts) and (arrows or cones):
            self._write_contacts(w, view, contacts, arrows, cones, size)
        self._write_legend(w, view, low, high)
        title = self._title(time, requested, len(ids), contacts, mode, arrows)
        w(f'<text x="6" y="-7" font-size="12" fill="{INK}">{_escape(title)}</text>\n')
        w("</svg>\n")

    def _view(self, centres, radii, geometry, size):
        """The bodies of the step and the static objects, the very long
        ones (Line) left out of the bounds."""
        points = [centres - radii[:, None], centres + radii[:, None]]
        if len(centres):
            extent = np.ptp(np.vstack(points), axis=0).max()
        else:
            extent = self._diagonal
        for cell in geometry.cells:
            p = geometry.points[np.asarray(cell)][:, :2]
            if np.ptp(p, axis=0).max() <= 100 * max(extent, size):
                points.append(p)
        points = np.vstack(points) if points else np.zeros((1, 2))
        low, high = points.min(axis=0), points.max(axis=0)
        diagonal = float(np.hypot(*(high - low)))
        return _View(low, high, max(size, 0.04 * diagonal))

    def _write_static(self, w, view, geometry, cell_ids, sides, size):
        w('<g id="static" stroke-linecap="butt">\n')
        diagonal = np.hypot(view.width, view.height) / view.scale
        step = max(HATCH_STEP * size, diagonal / HATCH_MAX)
        length = 0.9 * step
        rank = {}
        for cell, i in zip(geometry.cells, cell_ids):
            i = int(i)
            k = rank.get(i, 0)
            rank[i] = k + 1
            p = geometry.points[np.asarray(cell)][:, :2]
            if len(cell) > 2:
                w(
                    f'<polygon id="s{-i}_{k}" points="{view.points(p)}" '
                    f'fill="{STATIC_FILL}" stroke="{INK}" stroke-width="1.5"/>\n'
                )
                continue
            p1, p2 = p
            side = sides.get((i, k))
            if side is not None:
                # as the ground of the notes figures: a stroke going back
                # along the segment and into the material
                t = (p2 - p1) / np.linalg.norm(p2 - p1)
                d = (-0.2 * t + 0.25 * side) / 0.32 * length
                n = max(int(np.linalg.norm(p2 - p1) / step), 1)
                strokes = "".join(
                    f"M{view.xy(p1 + j / n * (p2 - p1))}l{view.delta(d)}"
                    for j in range(n + 1)
                )
                w(
                    f'<path id="h{-i}_{k}" d="{strokes}" stroke="{MUTED}" '
                    f'stroke-width="0.8" fill="none"/>\n'
                )
            w(
                f'<line id="s{-i}_{k}" x1="{view.x(p1)}" y1="{view.y(p1)}" '
                f'x2="{view.x(p2)}" y2="{view.y(p2)}" stroke="{INK}" '
                f'stroke-width="1.5"/>\n'
            )
        w("</g>\n")

    def _write_bodies(self, w, view, mode, ids, centres, angles, data, colours):
        blocks, scales = data["shape_index"], data["glyph_scale"]
        if mode == "bins":
            self._write_bins(w, view, centres, angles, blocks, scales, colours)
            return
        if mode == "use":
            # the shapes once, in the body frame
            w("<defs>\n")
            defined = {}
            for block, scale in sorted(set(zip(blocks.tolist(), scales.tolist()))):
                name = f"d{len(defined)}"
                defined[(block, scale)] = name
                if block == self._unit:
                    w(f'<circle id="{name}" r="{view.length(scale)}"/>\n')
                else:
                    geometry = self.shapes.blocks[block][1]
                    path = "".join(
                        "M"
                        + "L".join(
                            view.delta(q)
                            for q in scale * geometry.points[np.asarray(c)][:, :2]
                        )
                        + ("Z" if len(c) > 2 else "")
                        for c in geometry.cells
                    )
                    w(f'<path id="{name}" d="{path}"/>\n')
            w("</defs>\n")
        w(f'<g id="bodies" stroke="{INK}" stroke-width="0.8" fill-opacity="0.85">\n')
        for i, c, a, block, scale, k in zip(
            ids, centres, angles, blocks, scales, colours
        ):
            fill = RAMP[k]
            if mode == "use":
                rotation = (
                    f" rotate({-np.degrees(a):.3f})" if block != self._unit else ""
                )
                w(
                    f'<use id="b{i}" xlink:href="#{defined[(block, scale)]}" '
                    f'transform="translate({view.xy(c)}){rotation}" fill="{fill}"/>\n'
                )
            elif block == self._unit:
                tick = c + _rotate(a, np.array((scale, 0.0)))
                w(
                    f'<g id="b{i}"><circle cx="{view.x(c)}" cy="{view.y(c)}" '
                    f'r="{view.length(scale)}" fill="{fill}"/>'
                    f'<line x1="{view.x(c)}" y1="{view.y(c)}" x2="{view.x(tick)}" '
                    f'y2="{view.y(tick)}"/></g>\n'
                )
            else:
                geometry = self.shapes.blocks[block][1]
                w(f'<g id="b{i}" fill="{fill}">')
                for cell in geometry.cells:
                    p = c + _rotate(a, scale * geometry.points[np.asarray(cell)][:, :2])
                    if len(cell) > 2:
                        w(f'<polygon points="{view.points(p)}"/>')
                    else:
                        w(f'<polyline points="{view.points(p)}" fill="none"/>')
                w("</g>\n")
            if self.labels:
                w(
                    f'<text x="{view.x(c)}" y="{view.y(c)}" font-size="9" '
                    f'text-anchor="middle" fill="{INK}" stroke="none">{i}</text>\n'
                )
        w("</g>\n")

    def _write_bins(self, w, view, centres, angles, blocks, scales, colours):
        """The bodies merged in one path per colour, no outline."""
        w('<g id="bodies" stroke="none">\n')
        for k, fill in enumerate(RAMP):
            selected = np.flatnonzero(colours == k)
            if len(selected) == 0:
                continue
            w(f'<path fill="{fill}" d="')
            for j in selected:
                c, block, scale = centres[j], blocks[j], scales[j]
                if block == self._unit:
                    r = view.length(scale)
                    w(
                        f"M{view.xy(c - (scale, 0))}a{r},{r} 0 1 0 {2 * float(r):.2f},0"
                        f"a{r},{r} 0 1 0 {-2 * float(r):.2f},0"
                    )
                    continue
                geometry = self.shapes.blocks[block][1]
                for cell in geometry.cells:
                    if len(cell) > 2:
                        p = c + _rotate(
                            angles[j], scale * geometry.points[np.asarray(cell)][:, :2]
                        )
                        w(f"M{'L'.join(view.xy(q) for q in p)}Z")
            w('"/>\n')
        w("</g>\n")

    def _write_contacts(self, w, view, contacts, arrows, cones, size):
        """By (ds1, ds2, interaction): the Coulomb cone about the normal at the
        point of ds1, the reaction impulse on ds1 scaled to the largest."""
        order = np.lexsort(
            (
                contacts[:, CF["interaction"]],
                contacts[:, CF["ds2"]],
                contacts[:, CF["ds1"]],
            )
        )
        contacts = contacts[order]
        reactions = contacts[:, CF["reaction"]][:, :2]
        largest = float(np.linalg.norm(reactions, axis=1).max())
        length = 1.2 * size
        w('<g id="contacts" fill="none">\n')
        for row, r in zip(contacts, reactions):
            p = row[CF["point_b"]][:2]
            n = row[CF["normal"]][:2]
            if cones:
                half = np.arctan(row[CF["mu"]])
                e1 = p + length * _rotate(half, n)
                e2 = p + length * _rotate(-half, n)
                w(
                    f'<path d="M{view.xy(e1)}L{view.xy(p)}L{view.xy(e2)}" '
                    f'stroke="{CONE}" stroke-width="1" fill="{CONE}" fill-opacity="0.15"/>\n'
                )
            if arrows and largest > 0:
                tip = p + r / largest * length
                w(
                    f'<line x1="{view.x(p)}" y1="{view.y(p)}" x2="{view.x(tip)}" '
                    f'y2="{view.y(tip)}" stroke="{FORCE}" stroke-width="1.5"/>\n'
                )
        w("</g>\n")
        self._largest = largest

    def _write_legend(self, w, view, low, high):
        x, y = view.width - 80 - 18 * len(RAMP), view.height - 20
        w(f'<g id="legend" font-size="11" fill="{INK}">')
        w(
            f'<text x="{x - 4:.0f}" y="{y + 11:.0f}" text-anchor="end">|v| {low:.3g}</text>'
        )
        for k, colour in enumerate(RAMP):
            w(
                f'<rect x="{x + 18 * k:.0f}" y="{y:.0f}" width="18" height="12" fill="{colour}"/>'
            )
        w(
            f'<text x="{x + 18 * len(RAMP) + 4:.0f}" y="{y + 11:.0f}">'
            f"{high:.3g} m/s</text></g>\n"
        )

    def _title(self, time, requested, bodies, contacts, mode, arrows):
        title = f"t = {time:.6g} s"
        if requested is not None and requested != time:
            title += f" (requested {requested:.6g})"
        if self.scheme:
            title += f" · {self.scheme}"
        title += f" · {bodies} bodies · {len(contacts)} contacts"
        if len(contacts):
            mus = np.unique(np.round(contacts[:, CF["mu"]], 6))
            title += " · mu = " + ", ".join(f"{m:.3g}" for m in mus)
            if arrows and self._largest > 0:
                title += f" · arrows: reaction impulse, largest {self._largest:.3g} N s"
        if mode == "bins":
            title += " · bodies merged by colour"
        return title


# -- during a run ---------------------------------------------------------------


class RunPictures:
    """The pictures of the svg_output run option.

    requests: times (float) and simulation step numbers (int, the time k h).
    A request is drawn at the step k with t_k <= request < t_k + h, an output
    step or not. The votes for the sides are those of the output steps and of
    the drawn steps, as convert counts those of the file.
    """

    def __init__(self, io_filename, shapes, requests, h, directory=None, **options):
        self.io_filename = io_filename
        self.directory = directory
        self.h = h
        times = [
            k * h if isinstance(k, (int, np.integer)) else float(k) for k in requests
        ]
        self.pending = sorted(times)
        self.picture = Picture(shapes, **options)
        self._voted = None
        self.names = []

    def due(self, time):
        """The requests of the step at time."""
        eps = 1e-9 * self.h
        due = [t for t in self.pending if time - eps <= t < time + self.h - eps]
        return due

    def output_step(self, time, bodies, velocities, static_rows, contacts):
        """Count the votes of an output step; the first one gives the bodies of
        the segments without contact."""
        self.picture.first_bodies(bodies)
        if self._voted != time:
            self.picture.vote(time, static_rows, contacts)
            self._voted = time

    def draw(self, time, bodies, velocities, static_rows, contacts):
        """The pictures due at time, from the rows of the step (as_table)."""
        for requested in self.due(time):
            self.output_step(time, bodies, velocities, static_rows, contacts)
            name = svg_filename(self.io_filename, time, self.directory)
            with open(name, "w") as out:
                self.picture.write(
                    out, time, bodies, velocities, static_rows, contacts, requested
                )
            if abs(requested - time) > 1e-9 * self.h:
                warn(f"time {requested:.6g} requested, the step at {time:.6g} drawn")
            self.pending.remove(requested)
            self.names.append(name)

    def finish(self):
        for requested in self.pending:
            warn(f"time {requested:.6g} requested, not reached by the run: not drawn")


# -- from the HDF5 file -------------------------------------------------------


class FileSteps:
    """The output steps of a siconos mechanics HDF5 file, by its step
    index."""

    def __init__(self, data):
        if "step_index" not in data:
            raise ValueError("the file has no step index (data/step_index)")
        self.data = data
        index = data["step_index"]
        self.index = {name: index[name][:] for name in index}
        self.times = self.index["dynamic"][:, 0].astype(np.float64)
        self.contact_columns = data["cf"].shape[1] - 1 if "cf" in data else 25

    def resolve(self, time=None, frame=None):
        """The output step drawn for a time (the step at or before it, never
        interpolated) or a frame (an output step index, negative from the
        end); the last one if neither."""
        n = len(self.times)
        if frame is not None:
            if not -n <= frame < n:
                raise ValueError(f"frame {frame} not in 0..{n - 1}")
            return frame % n
        if time is None:
            return n - 1
        k = int(np.searchsorted(self.times, time, side="right")) - 1
        if k < 0:
            raise ValueError(
                f"time {time} before the first output step, {self.times[0]}"
            )
        if self.times[k] != time:
            warn(
                f"time {time:.6g} requested, the output step at {self.times[k]:.6g} drawn"
            )
        return k

    def rows(self, table, k, columns):
        """The rows of an output step without the time column."""
        if table not in self.index or table not in self.data:
            return np.zeros((0, columns))
        _, first, count = self.index[table][k]
        first, count = int(first), int(count)
        rows = self.data[table][first : first + count, 1:]
        return np.asarray(rows, dtype=np.float64).reshape(-1, columns)

    def step(self, k):
        """time, bodies, velocities, static rows, contacts of output step k."""
        return (
            float(self.times[k]),
            self.rows("dynamic", k, 8),
            self.rows("velocities", k, 7),
            self.rows("static", k, 8),
            self.rows("cf", k, self.contact_columns),
        )


def convert(io_filename, times=(), frames=(), directory=None, output=None, **options):
    """Write the pictures of the output steps at times (at or before each)
    and frames (output step indices) of a siconos mechanics HDF5 file, the
    last step if none. options: those of Picture. Returns the file names.
    """
    import h5py

    names = []
    with h5py.File(io_filename, "r") as f:
        dimension = int(f.attrs.get("dimension", 3))
        if dimension != 2:
            raise ValueError(f"{io_filename}: SVG pictures are for 2D files")
        data = f["data"]
        steps = FileSteps(data)
        options.setdefault("scheme", _file_scheme(data))
        picture = Picture(vtkhdf.ShapeTable(data, dimension), **options)
        requests = [(steps.resolve(time=t), t) for t in times]
        requests += [(steps.resolve(frame=k), None) for k in frames]
        if not requests:
            requests = [(steps.resolve(), None)]
        if output is not None and len(requests) > 1:
            raise ValueError("one output file name for several pictures")
        # the votes of the steps up to each drawn one, in order
        voted = -1
        first = steps.step(0)
        picture.first_bodies(first[1])
        for k, requested in sorted(requests, key=lambda r: r[0]):
            for j in range(voted + 1, k + 1):
                time, _, _, static_rows, contacts = steps.step(j)
                picture.vote(time, static_rows, contacts)
            voted = max(voted, k)
            time, bodies, velocities, static_rows, contacts = steps.step(k)
            name = output or svg_filename(io_filename, time, directory)
            with open(name, "w") as out:
                picture.write(
                    out, time, bodies, velocities, static_rows, contacts, requested
                )
            names.append(name)
    return names
