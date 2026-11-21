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
from collections import namedtuple

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
# the title: lines of at most TITLE_CHARACTERS characters, cut at the " · "
TITLE_LINE = 16
TITLE_CHARACTERS = 110

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
# a contact at an end point of a segment does not vote, nor one with a normal
# closer to the segment direction than this cosine
END_POINT_COSINE = 0.1
# two times closer than this, relative, are the same step (the times of the
# tables are float32: 6e-8 relative)
SAME_TIME = 1e-6


def _same_time(a, b):
    return abs(a - b) <= SAME_TIME * max(1.0, abs(a))


def scheme_name(osi):
    """The name of an integrator (a class or an object) for the title, as
    the run options of the file have it."""
    if osi is None or isinstance(osi, str):
        return osi
    return getattr(osi, "__name__", type(osi).__name__)


def _file_run_options(data):
    """The run options of a file, {} if not written."""
    import json

    options = data.get("siconos_mechanics_run_options")
    if options is None or "options" not in options.attrs:
        return {}
    return json.loads(options.attrs["options"])


def _units(length_unit):
    """The units of a velocity and of an impulse: SI, or lengths in u =
    length_unit m (a run with a gravity_scale)."""
    if length_unit in (None, 1, 1.0):
        return "m/s", "N s"
    return "u/s", "kg u/s"


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
            distance, t = _distance_to_segment(points, p1, p2)
            # a contact at an end point (a corner of a polyline, against the
            # face of a polygon) has the normal of the other body: no vote
            end = self.tolerance / np.linalg.norm(p2 - p1)
            on = (distance < self.tolerance) & (t > end) & (t < 1 - end)
            dots = normals[on] @ _perp(p1, p2)
            dots = dots[np.abs(dots) > END_POINT_COSINE]
            if len(dots):
                count = self.votes.setdefault(key, [0, 0])
                count[0] += int(np.sum(dots > 0))
                count[1] += int(np.sum(dots < 0))

    def _raw(self, key, p1, p2, centre_sets, size):
        """(votes on the side of perp, on the other, where they come from)."""
        plus, minus = self.votes.get(key, (0, 0))
        if plus + minus:
            return plus, minus, "contacts"
        perp = _perp(p1, p2)
        for centres in centre_sets:
            distance, t = _distance_to_segment(centres, p1, p2)
            near = (distance < 2 * size) & (t > 0) & (t < 1)
            side = (centres[near] - p1) @ perp
            plus, minus = int(np.sum(side > 0)), int(np.sum(side < 0))
            if plus + minus:
                return plus, minus, "bodies near it"
        return 0, 0, None

    def _chains(self, segments):
        """The segments joined end to end, two at a joint (not three or more:
        a wall on a floor is a junction), as lists of (index, aligned):
        aligned when the segment runs in the direction of its chain."""
        ends = {}
        for k, (_, p1, p2) in enumerate(segments):
            for end, p in enumerate((p1, p2)):
                q = tuple(np.round(np.asarray(p) / self.tolerance).astype(np.int64))
                ends.setdefault(q, []).append((k, end))
        neighbours = {k: [] for k in range(len(segments))}
        for joint in ends.values():
            if len(joint) == 2 and joint[0][0] != joint[1][0]:
                (a, ea), (b, eb) = joint
                neighbours[a].append((b, ea != eb))
                neighbours[b].append((a, ea != eb))
        chains, seen = [], set()
        for k in range(len(segments)):
            if k in seen:
                continue
            chain, stack = [], [(k, True)]
            seen.add(k)
            while stack:
                a, aligned = stack.pop()
                chain.append((a, aligned))
                for b, same in neighbours[a]:
                    if b not in seen:
                        seen.add(b)
                        stack.append((b, aligned if same else not aligned))
            chains.append(chain)
        return chains

    def material(self, segments, centre_sets, size):
        """key -> unit normal pointing into the material, or None (not
        hatched), with a warning for the segments not hatched.

        The votes of a segment are its contacts, else the bodies near it
        (centre_sets: the body centres to try, in order). A chain of
        segments with votes on one side only takes that side, its segments
        without vote included (a terrain the flow has not reached); else
        each segment by its own votes.
        """
        raw = [self._raw(key, p1, p2, centre_sets, size) for key, p1, p2 in segments]
        sides, none, both = {}, [], []
        for chain in self._chains(segments):
            # the votes in the direction of the chain
            plus = sum(raw[k][0] if a else raw[k][1] for k, a in chain)
            minus = sum(raw[k][1] if a else raw[k][0] for k, a in chain)
            one_side = len(chain) > 1 and (plus > 0) != (minus > 0)
            for k, aligned in chain:
                key, p1, p2 = segments[k]
                perp = _perp(p1, p2)
                if key[0] in self.forced:
                    sides[key] = self.forced[key[0]] * perp
                    continue
                if one_side:
                    bodies = 1 if (plus > 0) == aligned else -1
                    sides[key] = -bodies * perp
                    continue
                p, m, source = raw[k]
                if p + m == 0:
                    sides[key] = None
                    none.append((key, p1, p2))
                elif p and m:
                    sides[key] = None
                    both.append((key, p1, p2, source, p, m))
                else:
                    # the material is on the other side of the bodies
                    sides[key] = -perp if p else perp
        for key, p1, p2, source, p, m in both:
            warn(
                f"static {key[0]} segment {_point(p1)} {_point(p2)}: bodies on "
                f"both sides ({source}: {p} and {m}), not hatched"
            )
        if len(none) <= 5:
            for key, p1, p2 in none:
                warn(
                    f"static {key[0]} segment {_point(p1)} {_point(p2)}: no "
                    f"contact nor body near it, not hatched"
                )
        else:
            ids = ", ".join(str(key[0]) for key, _, _ in none[:10])
            warn(
                f"{len(none)} segments with no contact nor body near them, not "
                f"hatched: static {ids}{', ...' if len(none) > 10 else ''}"
            )
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


FEM_FILL = "#cfe0ee"
# the bodies when they are not coloured
BODY = "#c9d6e3"
# the colourings of the bodies: name, legend, unit
COLOURINGS = {"speed": ("|v|", "m/s"), "omega": ("|ω|", "rad/s"), "none": None}
# the force chains: widths in px of the strongest and of the weakest, in
# CHAIN_BINS classes
CHAIN_WIDTHS = (0.5, 4.0)
CHAIN_BINS = 8
# the arrow scales: proportional to the impulse, a fixed length (the
# direction only), or logarithmic over LOG_DECADES decades below the largest
SCALES = ("proportional", "fixed", "log")
LOG_DECADES = 3
# a contact is at the limit of its cone if |r_t| >= (1 - SATURATION) mu r_n
SATURATION = 1e-3
# the animations: the duration of a frame when one frame per output step
FRAME_DURATION = 0.1


class Step(
    namedtuple("Step", "time bodies velocities static contacts fem", defaults=(None,))
):
    """The rows of an output step, as the HDF5 tables have them without their
    time column: bodies (N, 8) id, x, y, z, q0..q3; velocities (N, 7) id, v,
    w; static (M, 8) as bodies; contacts, the cf rows; fem, {ds id: the
    displacements of the nodes in spatial order} (empty: undeformed)."""


class _Frame:
    """What a step draws: the bodies by id, their pose and colour value, the
    static objects, the contacts in their order."""

    def __init__(self, picture, step):
        bodies = step.bodies[np.argsort(step.bodies[:, 0], kind="stable")]
        velocities = step.velocities[np.argsort(step.velocities[:, 0], kind="stable")]
        self.step = step
        self.ids = bodies[:, 0].astype(np.int64)
        self.centres = bodies[:, 1:3]
        self.angles = 2 * np.arctan2(bodies[:, 7], bodies[:, 4])
        self.data = picture.shapes.point_data(self.ids)
        if picture.colour == "omega":
            self.values = np.abs(velocities[:, 6])
        else:
            self.values = np.linalg.norm(velocities[:, 1:3], axis=1)
        self.geometry, self.cell_ids, self.segments = picture.segments(
            step.static, step.time
        )
        contacts = step.contacts
        order = np.lexsort(
            (
                contacts[:, CF["interaction"]],
                contacts[:, CF["ds2"]],
                contacts[:, CF["ds1"]],
            )
        )
        self.contacts = contacts[order]
        self.fem = picture.fem_points(step.fem)

    def extent_points(self, picture):
        """The points the view must hold."""
        radii = self.data["disk_radius"][:, None]
        points = [self.centres - radii, self.centres + radii]
        extent = (
            np.ptp(np.vstack(points), axis=0).max()
            if len(self.centres)
            else picture._diagonal
        )
        for cell in self.geometry.cells:
            p = self.geometry.points[np.asarray(cell)][:, :2]
            # the very long static objects (Line) left out
            if np.ptp(p, axis=0).max() <= 100 * max(
                extent, picture.shapes.typical_size
            ):
                points.append(p)
        if self.fem is not None:
            points.append(self.fem)
        return np.vstack(points) if points else np.zeros((1, 2))


class Picture:
    """The pictures of a simulation: its shapes (a vtkhdf.ShapeTable), its
    FEM meshes (a vtkhdf.FemMeshes, or None), the votes for the sides of its
    segments, the drawing options.

    mode: "auto", "circle", "use" or "bins" (see CIRCLE_LIMIT).
    arrows, chains, cones: True, False or None (by the number of contacts):
    the reaction impulse on ds1 at its point, the force chains (a line
    between the bodies, its width by r_n), the Coulomb cones.
    scale: the length of the arrows, "proportional" to the impulse, "fixed"
    (the direction only) or "log"; the edges of a cone are as long as the
    arrow, so the tip of the arrow is on an edge when the contact slides.
    colour: "speed", "omega" or "none".
    material_sides: {static id: +1 | -1} the material side of the
    segments of a static object, on the side of perp(p2 - p1) (+1) or
    the other, instead of the inference.
    scheme: the name of the integrator, in the title.
    """

    def __init__(
        self,
        shapes,
        fem=None,
        mode="auto",
        arrows=None,
        chains=None,
        cones=None,
        scale="proportional",
        colour="speed",
        material_sides=None,
        labels=False,
        scheme=None,
        length_unit=1.0,
        view=None,
    ):
        if shapes.dimension != 2:
            raise ValueError("SVG pictures are for 2D simulations")
        if scale not in SCALES:
            raise ValueError(f"scale {scale!r}: one of {', '.join(SCALES)}")
        if colour not in COLOURINGS:
            raise ValueError(f"colour {colour!r}: one of {', '.join(COLOURINGS)}")
        if mode not in ("auto", "circle", "use", "bins"):
            raise ValueError(f"mode {mode!r}: auto, circle, use or bins")
        self.shapes = shapes
        self.fem = fem if fem else None
        self.mode = mode
        self.arrows = arrows
        self.chains = chains
        self.cones = cones
        self.scale = scale
        self.colour = colour
        self.labels = labels
        self.scheme = scheme
        self.length_unit = 1.0 if length_unit is None else float(length_unit)
        # a crop: x0, x1, y0, y1
        self.view = None
        if view is not None:
            x0, x1, y0, y1 = map(float, view)
            self.view = (np.array((x0, y0)), np.array((x1, y1)))
        bounds = shapes.bounds or [0, 1, 0, 1, 0, 0]
        self._diagonal = float(np.hypot(bounds[1] - bounds[0], bounds[3] - bounds[2]))
        self.sides = Sides(ON_SEGMENT * self._diagonal, material_sides)
        self._unit = shapes._block_of_key.get(("unit",))
        self._first_centres = None

    # -- the steps

    def segments(self, static_rows, time):
        geometry, ids = self.shapes.static_geometry(static_rows, time)
        return geometry, ids, _segments(geometry, ids)

    def vote(self, step):
        """Count the contacts of an output step in the votes."""
        if len(step.contacts):
            _, _, segments = self.segments(step.static, step.time)
            self.sides.add_contacts(segments, step.contacts)

    def first_bodies(self, bodies):
        """The bodies of the first output step, for the segments without
        contact."""
        if self._first_centres is None:
            self._first_centres = np.asarray(bodies)[:, 1:3].copy()

    def fem_points(self, displacements):
        """The deformed nodes of the FEM meshes, None without mesh."""
        if self.fem is None:
            return None
        points, _, _ = self.fem.step_data(displacements or {}, {}, {})
        return points[:, :2]

    def _sides(self, frame, centres=None):
        sets = [frame.centres if centres is None else centres]
        if self._first_centres is not None:
            sets.append(self._first_centres)
        return self.sides.material(frame.segments, sets, self.shapes.typical_size)

    def _mode(self, n):
        if self.mode != "auto":
            return self.mode
        return "circle" if n <= CIRCLE_LIMIT else "use" if n <= USE_LIMIT else "bins"

    def _contact_layers(self, n):
        arrows = n <= ARROW_LIMIT if self.arrows is None else self.arrows
        chains = n > ARROW_LIMIT if self.chains is None else self.chains
        cones = n <= CONE_LIMIT if self.cones is None else self.cones
        return arrows, chains, cones

    def bounds(self, frames):
        """(low, high) of the view of frames."""
        points = np.vstack([f.extent_points(self) for f in frames])
        return points.min(axis=0), points.max(axis=0)

    def _view(self, bounds):
        low, high = bounds
        diagonal = float(np.hypot(*(high - low)))
        return _View(low, high, max(self.shapes.typical_size, 0.04 * diagonal))

    # -- a still picture

    def write(self, out, step, requested=None, bounds=None):
        """Write the picture of an output step in the text file out. The votes
        of the steps up to this one must have been counted (vote). bounds:
        those of an animation (bounds), for its frame; else the step's."""
        frame = _Frame(self, step)
        view = self._view(self.view or bounds or self.bounds([frame]))
        low, high = self._range([frame])
        mode = self._mode(len(frame.ids))
        arrows, chains, cones = self._contact_layers(len(frame.contacts))
        largest = self._largest(frame)
        title = self._title(step.time, requested, frame, mode, arrows, chains, largest)
        w = out.write
        self._header(w, view, title)
        self._write_static(w, view, frame, self._sides(frame))
        self._write_fem(w, view, frame.fem)
        self._write_bodies(w, view, mode, frame, low, high)
        self._write_contacts(w, view, frame, arrows, chains, cones)
        self._write_legend(w, view, low, high)
        w("</svg>\n")

    def _largest(self, frame):
        """The largest reaction impulse of a frame."""
        if len(frame.contacts) == 0:
            return 0.0
        return self._arrow_lengths(frame.contacts[:, CF["reaction"]][:, :2])[1]

    def _header(self, w, view, title):
        """The svg element, the background, the arrow head and the title,
        above the view."""
        lines = [""]
        for part in title.split(" · "):
            if lines[-1] and len(lines[-1]) + 3 + len(part) > TITLE_CHARACTERS:
                lines.append(part)
            else:
                lines[-1] += (" · " if lines[-1] else "") + part
        top = TITLE_LINE * len(lines) + 6
        height = view.height + top
        w('<?xml version="1.0" encoding="UTF-8"?>\n')
        w(
            f'<svg xmlns="http://www.w3.org/2000/svg" '
            f'xmlns:xlink="http://www.w3.org/1999/xlink" '
            f'width="{view.width:.0f}" height="{height:.0f}" '
            f'viewBox="0 {-top} {view.width:.0f} {height:.0f}" '
            f'font-family="sans-serif">\n'
        )
        w(f'<rect x="0" y="{-top}" width="100%" height="100%" fill="{BACKGROUND}"/>\n')
        w(f'<text id="title" font-size="12" fill="{INK}">')
        for k, line in enumerate(lines):
            y = -top + TITLE_LINE * (k + 1)
            w(f'<tspan x="6" y="{y}">{_escape(line)}</tspan>')
        w("</text>\n")
        w(
            f'<defs><marker id="head" viewBox="0 0 10 10" refX="9" refY="5" '
            f'markerWidth="5" markerHeight="5" orient="auto-start-reverse">'
            f'<path d="M0,0L10,5L0,10z" fill="{FORCE}"/></marker></defs>\n'
        )
        if self.view is not None:
            # a crop: the scene stays below the title, closed by the legend
            w(
                f'<clipPath id="crop"><rect x="0" y="0" width="{view.width:.0f}" '
                f'height="{view.height:.0f}"/></clipPath>\n'
                f'<g clip-path="url(#crop)">\n'
            )

    def _range(self, frames):
        values = [f.values for f in frames if len(f.values)]
        if self.colour == "none" or not values:
            return 0.0, 0.0
        values = np.concatenate(values)
        return float(values.min()), float(values.max())

    def _fills(self, values, low, high):
        if self.colour == "none":
            return [BODY] * len(values)
        return [RAMP[k] for k in _ramp_index(values, low, high)]

    def _write_static(self, w, view, frame, sides, animate=None):
        """One group per static object, id s<-id>: its polygons filled, its
        segments with their hatch. animate: {static id: the <animate> of its
        visibility}."""
        geometry, cell_ids = frame.geometry, frame.cell_ids
        diagonal = np.hypot(view.width, view.height) / view.scale
        step = max(HATCH_STEP * self.shapes.typical_size, diagonal / HATCH_MAX)
        length = 0.9 * step
        w('<g id="static" stroke-linecap="butt">\n')
        rank, current = {}, None
        for cell, i in zip(geometry.cells, cell_ids):
            i = int(i)
            if i != current:
                if current is not None:
                    w("</g>\n")
                w(f'<g id="s{-i}">')
                if animate and i in animate:
                    w(animate[i])
                current = i
            k = rank.get(i, 0)
            rank[i] = k + 1
            p = geometry.points[np.asarray(cell)][:, :2]
            if len(cell) > 2:
                w(
                    f'<polygon id="p{-i}_{k}" points="{view.points(p)}" '
                    f'fill="{STATIC_FILL}" stroke="{INK}" stroke-width="1.5"/>'
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
                    f'stroke-width="0.8" fill="none"/>'
                )
            w(
                f'<line id="l{-i}_{k}" x1="{view.x(p1)}" y1="{view.y(p1)}" '
                f'x2="{view.x(p2)}" y2="{view.y(p2)}" stroke="{INK}" '
                f'stroke-width="1.5"/>'
            )
        if current is not None:
            w("</g>\n")
        w("</g>\n")

    def _write_fem(self, w, view, points, animate=None):
        """The deformed meshes, one path per body. animate: the paths of the
        other frames, {ds id: [d, ...]}, and the timing of their <animate>."""
        if points is None:
            return
        w(f'<g id="fem" fill="{FEM_FILL}" stroke="{INK}" stroke-width="0.4">\n')
        for ds_id in self.fem.ids:
            cells = self.fem.cell_ranges[ds_id]
            triangles = self.fem.triangles[cells]
            d = "".join(
                f"M{view.xy(points[a])}L{view.xy(points[b])}L{view.xy(points[c])}Z"
                for a, b, c in triangles
            )
            w(f'<path id="f{ds_id}" d="{d}">')
            if animate:
                paths, timing = animate
                w(
                    f'<animate attributeName="d" values="{";".join(paths[ds_id])}"{timing}/>'
                )
            w("</path>\n")
        w("</g>\n")

    def _shape(self, w, view, block, scale, centre=None, angle=0.0):
        """A body shape, at its pose, or in its frame (centre None)."""
        if block == self._unit:
            if centre is None:
                w(
                    f'<circle r="{view.length(scale)}"/>'
                    f'<line x2="{view.length(scale)}" y2="0"/>'
                )
            else:
                tick = centre + _rotate(angle, np.array((scale, 0.0)))
                w(
                    f'<circle cx="{view.x(centre)}" cy="{view.y(centre)}" '
                    f'r="{view.length(scale)}"/>'
                    f'<line x1="{view.x(centre)}" y1="{view.y(centre)}" '
                    f'x2="{view.x(tick)}" y2="{view.y(tick)}"/>'
                )
            return
        geometry = self.shapes.blocks[block][1]
        for cell in geometry.cells:
            q = scale * geometry.points[np.asarray(cell)][:, :2]
            if centre is None:
                points = " ".join(view.delta(p) for p in q)
            else:
                points = view.points(centre + _rotate(angle, q))
            if len(cell) > 2:
                w(f'<polygon points="{points}"/>')
            else:
                w(f'<polyline points="{points}" fill="none"/>')

    def _write_bodies(self, w, view, mode, frame, low, high, opacity=0.85, group=True):
        ids, centres, angles = frame.ids, frame.centres, frame.angles
        blocks, scales = frame.data["shape_index"], frame.data["glyph_scale"]
        fills = self._fills(frame.values, low, high)
        if mode == "bins":
            self._write_bins(w, view, frame, fills)
            return
        if mode == "use":
            # the shapes once, in the body frame, without the rotation tick
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
        if group:
            w(
                f'<g id="bodies" stroke="{INK}" stroke-width="0.8" '
                f'fill-opacity="{opacity}">\n'
            )
        for i, c, a, block, scale, fill in zip(
            ids, centres, angles, blocks, scales, fills
        ):
            if mode == "use":
                rotation = (
                    f" rotate({-np.degrees(a):.3f})" if block != self._unit else ""
                )
                w(
                    f'<use id="b{i}" xlink:href="#{defined[(block, scale)]}" '
                    f'transform="translate({view.xy(c)}){rotation}" fill="{fill}"/>\n'
                )
            else:
                w(f'<g id="b{i}" fill="{fill}">')
                self._shape(w, view, block, scale, c, a)
                w("</g>\n")
            if self.labels:
                w(
                    f'<text x="{view.x(c)}" y="{view.y(c)}" font-size="9" '
                    f'text-anchor="middle" fill="{INK}" stroke="none">{i}</text>\n'
                )
        if group:
            w("</g>\n")

    def _write_bins(self, w, view, frame, fills):
        """The bodies merged in one path per colour, no outline."""
        centres, angles = frame.centres, frame.angles
        blocks, scales = frame.data["shape_index"], frame.data["glyph_scale"]
        fills = np.asarray(fills)
        w('<g id="bodies" stroke="none">\n')
        for fill in dict.fromkeys(fills.tolist()):
            w(f'<path fill="{fill}" d="')
            for j in np.flatnonzero(fills == fill):
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

    def _arrow_lengths(self, reactions):
        """The lengths of the arrows of reactions, and the largest impulse."""
        norms = np.linalg.norm(reactions, axis=1)
        largest = float(norms.max()) if len(norms) else 0.0
        length = 1.2 * self.shapes.typical_size
        if largest <= 0:
            return np.zeros(len(norms)), largest
        if self.scale == "fixed":
            lengths = np.where(norms > 0, length, 0.0)
        elif self.scale == "log":
            floor = largest * 10.0**-LOG_DECADES
            lengths = (
                length * np.log10(1 + norms / floor) / np.log10(1 + 10**LOG_DECADES)
            )
        else:
            lengths = norms / largest * length
        return lengths, largest

    def _write_contacts(
        self,
        w,
        view,
        frame,
        arrows,
        chains,
        cones,
        group_id="contacts",
        extra="",
        follow=None,
    ):
        """The force chains, the cones and the arrows of a frame; returns the
        largest impulse. follow(ds1): the opening and closing markup of a
        group of the contacts of the body ds1 (an animation moves them with
        the body), None for one group."""
        contacts = frame.contacts
        if len(contacts) == 0 or not (arrows or chains or cones):
            return 0.0
        reactions = contacts[:, CF["reaction"]][:, :2]
        lengths, largest = self._arrow_lengths(reactions)
        w(f'<g id="{group_id}" fill="none"{extra}>\n')
        if follow is None:
            groups = [(None, np.arange(len(contacts)))]
        else:
            # the contacts are sorted by ds1 first
            ds1 = contacts[:, CF["ds1"]].astype(np.int64)
            bodies, starts = np.unique(ds1, return_index=True)
            ends = np.append(starts[1:], len(contacts))
            groups = [
                (int(b), np.arange(s, e)) for b, s, e in zip(bodies, starts, ends)
            ]
        size = self.shapes.typical_size
        for body, rows in groups:
            opening, closing = ("", "") if follow is None else follow(body)
            w(opening)
            if chains:
                self._write_chains(w, view, frame, rows)
            for j in rows:
                row, r, length = contacts[j], reactions[j], lengths[j]
                p = row[CF["point_b"]][:2]
                n = row[CF["normal"]][:2]
                if cones:
                    half = np.arctan(row[CF["mu"]])
                    edge = length if arrows and length > 0 else 0.6 * size
                    e1 = p + edge * _rotate(half, n)
                    e2 = p + edge * _rotate(-half, n)
                    w(
                        f'<path d="M{view.xy(e1)}L{view.xy(p)}L{view.xy(e2)}" '
                        f'stroke="{CONE}" stroke-width="1" fill="{CONE}" '
                        f'fill-opacity="0.15"/>\n'
                    )
                if arrows and length > 0:
                    tip = p + r / np.linalg.norm(r) * length
                    w(
                        f'<line x1="{view.x(p)}" y1="{view.y(p)}" x2="{view.x(tip)}" '
                        f'y2="{view.y(tip)}" stroke="{FORCE}" stroke-width="1.5" '
                        f'marker-end="url(#head)"/>\n'
                    )
            w(closing)
        w("</g>\n")
        return largest

    def _write_chains(self, w, view, frame, rows=None):
        """A line between the centres of the bodies of a contact (to the
        contact point for an obstacle), its width by r_n, one path per width
        class; the contacts of rows (all by default), the widths by the
        strongest of the frame."""
        contacts = frame.contacts
        rn = contacts[:, CF["reaction_local"]][:, 0]
        strongest = float(rn.max()) if len(rn) else 0.0
        if strongest <= 0:
            return
        selection = np.zeros(len(contacts), dtype=bool)
        selection[np.arange(len(contacts)) if rows is None else rows] = True
        ids = frame.ids

        def centre(body):
            k = np.searchsorted(ids, body)
            k = np.minimum(k, len(ids) - 1)
            return frame.centres[k], ids[k] == body

        c1, found1 = centre(contacts[:, CF["ds1"]].astype(np.int64))
        c2, found2 = centre(contacts[:, CF["ds2"]].astype(np.int64))
        obstacle = contacts[:, CF["ds1"]] == contacts[:, CF["ds2"]]
        c2 = np.where(obstacle[:, None], contacts[:, CF["point_a"]][:, :2], c2)
        keep = (rn > 0) & found1 & found2 & selection
        classes = np.rint(rn / strongest * (CHAIN_BINS - 1)).astype(np.int64)
        low, high = CHAIN_WIDTHS
        for k in range(CHAIN_BINS):
            selected = np.flatnonzero(keep & (classes == k))
            if len(selected) == 0:
                continue
            width = low + (high - low) * k / (CHAIN_BINS - 1)
            d = "".join(f"M{view.xy(c1[j])}L{view.xy(c2[j])}" for j in selected)
            w(
                f'<path d="{d}" stroke="{FORCE}" stroke-opacity="0.7" '
                f'stroke-width="{width:.2f}" stroke-linecap="round"/>\n'
            )

    def _write_legend(self, w, view, low, high):
        if self.view is not None:
            w("</g>\n")
        if self.colour == "none":
            return
        name, unit = COLOURINGS[self.colour]
        if unit == "m/s":
            unit = _units(self.length_unit)[0]
        x, y = view.width - 80 - 18 * len(RAMP), view.height - 20
        if high <= low:
            # one value (a single body, or bodies at rest): no ramp
            w(
                f'<g id="legend" font-size="11" fill="{INK}">'
                f'<rect x="{x:.0f}" y="{y:.0f}" width="18" height="12" '
                f'fill="{RAMP[0]}"/>'
                f'<text x="{x + 22:.0f}" y="{y + 11:.0f}">'
                f"{name} = {low:.3g} {unit}</text></g>\n"
            )
            return
        w(f'<g id="legend" font-size="11" fill="{INK}">')
        w(
            f'<text x="{x - 4:.0f}" y="{y + 11:.0f}" text-anchor="end">'
            f"{name} {low:.3g}</text>"
        )
        for k, colour in enumerate(RAMP):
            w(
                f'<rect x="{x + 18 * k:.0f}" y="{y:.0f}" width="18" height="12" '
                f'fill="{colour}"/>'
            )
        w(
            f'<text x="{x + 18 * len(RAMP) + 4:.0f}" y="{y + 11:.0f}">'
            f"{high:.3g} {unit}</text></g>\n"
        )

    def _contacts_title(self, contacts, arrows, chains, largest):
        title = f" · {len(contacts)} contacts"
        if len(contacts):
            mus = np.unique(np.round(contacts[:, CF["mu"]], 6))
            title += " · mu = " + ", ".join(f"{m:.3g}" for m in mus)
        if arrows and largest > 0:
            what = {
                "proportional": "reaction impulse",
                "fixed": "direction of the reaction impulse",
                "log": f"reaction impulse, log scale over {LOG_DECADES} decades",
            }[self.scale]
            title += (
                f" · arrows: {what}, largest {largest:.3g} "
                f"{_units(self.length_unit)[1]}"
            )
        if chains:
            title += " · chains: width by r_n"
        return title

    def _title(self, time, requested, frame, mode, arrows, chains, largest):
        title = f"t = {time:.6g} s"
        if requested is not None and not _same_time(requested, time):
            title += f" (requested {requested:.6g})"
        if self.scheme:
            title += f" · {self.scheme}"
        title += f" · {len(frame.ids)} bodies"
        if self.fem is not None:
            title += f" · {len(self.fem.ids)} FEM bodies"
        title += self._contacts_title(frame.contacts, arrows, chains, largest)
        title += self._view_title()
        if mode == "bins":
            title += " · bodies merged by colour"
        return title

    def _view_title(self):
        title = ""
        if self.length_unit != 1.0:
            title += f" · u = {self.length_unit:g} m"
        if self.view is not None:
            (x0, y0), (x1, y1) = self.view
            title += f" · view cropped to [{x0:g}, {x1:g}] x [{y0:g}, {y1:g}]"
        return title

    # -- several steps

    def write_strobe(self, out, steps, bounds=None):
        """The bodies of steps superimposed, the earlier ones fainter; the
        static objects and the contacts of the last step."""
        frames = [_Frame(self, s) for s in steps]
        last = frames[-1]
        view = self._view(self.view or bounds or self.bounds(frames))
        low, high = self._range(frames)
        # no <defs> per frame: one element per body, or the colour paths
        mode = "bins" if self._mode(len(last.ids)) == "bins" else "circle"
        arrows, chains, cones = self._contact_layers(len(last.contacts))
        n = len(frames)
        title = (
            f"t = {frames[0].step.time:.6g} to {last.step.time:.6g} s, "
            f"{n} steps superimposed"
        )
        if self.scheme:
            title += f" · {self.scheme}"
        title += f" · {len(last.ids)} bodies"
        title += self._contacts_title(
            last.contacts, arrows, chains, self._largest(last)
        )
        title += self._view_title()
        w = out.write
        self._header(w, view, title)
        self._write_static(w, view, last, self._sides(last))
        self._write_fem(w, view, last.fem)
        w(f'<g id="bodies" stroke="{INK}" stroke-width="0.8">\n')
        for k, frame in enumerate(frames):
            opacity = 0.15 + 0.85 * k / max(n - 1, 1)
            w(f'<g id="frame{k}" opacity="{opacity:.3f}">\n')
            self._write_bodies(w, view, mode, frame, low, high, group=False)
            w("</g>\n")
        w("</g>\n")
        self._write_contacts(w, view, last, arrows, chains, cones)
        self._write_legend(w, view, low, high)
        w("</svg>\n")

    def write_animation(
        self, out, steps, uniform=False, discrete=False, speed=1.0, contacts=None
    ):
        """An animated SVG (SMIL) of the steps: the bodies move, the static
        objects are shown while alive, the contacts of each step in turn. A
        viewer without SMIL shows the first step.

        uniform: one frame per step, FRAME_DURATION each, else the physical
        times of the steps, slowed down by speed (speed 1: real time).
        discrete: no interpolation between the steps.
        contacts: draw the contacts (None: below ARROW_LIMIT contacts in all
        the steps).
        The votes of all the steps must have been counted.
        """
        frames = [_Frame(self, s) for s in steps]
        n = len(frames)
        times = np.array([f.step.time for f in frames])
        view = self._view(self.view or self.bounds(frames))
        low, high = self._range(frames)
        if uniform or n < 2 or times[-1] == times[0]:
            key_times = np.linspace(0, 1, n) if n > 1 else np.zeros(1)
            duration = FRAME_DURATION * n
        else:
            key_times = (times - times[0]) / (times[-1] - times[0])
            duration = (times[-1] - times[0]) / speed
        mode = "discrete" if discrete else "linear"
        keys = ";".join(f"{t:.6f}" for t in key_times)
        timing = (
            f' keyTimes="{keys}" dur="{duration:.6g}s" calcMode="{mode}" '
            f'repeatCount="indefinite"'
        )
        discrete_timing = (
            f' keyTimes="{keys}" dur="{duration:.6g}s" calcMode="discrete" '
            f'repeatCount="indefinite"'
        )

        total = sum(len(f.contacts) for f in frames)
        most = max(len(f.contacts) for f in frames)
        draw = most <= ARROW_LIMIT if contacts is None else contacts
        playback = (
            "one frame per output step"
            if uniform
            else ("real time" if speed == 1 else f"{speed:.3g} x real time")
        )
        title = f"t = {times[0]:.6g} to {times[-1]:.6g} s, {n} output steps, {playback}"
        if self.scheme:
            title += f" · {self.scheme}"
        title += f" · {len(frames[0].ids)} bodies"
        if draw and total:
            title += " · contacts of each step"
        title += self._view_title()
        w = out.write
        self._header(w, view, title)

        # the static objects of all the steps, shown while alive; the sides
        # are those of all the steps
        static_ids = []
        alive = {}
        geometry_frame = None
        for k, frame in enumerate(frames):
            for i in np.unique(frame.cell_ids).tolist():
                if i not in alive:
                    alive[i] = np.zeros(n, dtype=bool)
                    static_ids.append(i)
                alive[i][k] = True
        if static_ids:
            # the frame with the most static objects gives the geometry
            geometry_frame = max(frames, key=lambda f: len(np.unique(f.cell_ids)))
            animate = {
                i: (
                    '<animate attributeName="visibility" values="'
                    + ";".join("visible" if a else "hidden" for a in alive[i])
                    + f'"{discrete_timing}/>'
                )
                for i in static_ids
                if not alive[i].all()
            }
            self._write_static(
                w, view, geometry_frame, self._sides(geometry_frame), animate
            )

        if self.fem is not None:
            paths = {
                ds_id: [
                    "".join(
                        f"M{view.xy(f.fem[a])}L{view.xy(f.fem[b])}L{view.xy(f.fem[c])}Z"
                        for a, b, c in self.fem.triangles[self.fem.cell_ranges[ds_id]]
                    )
                    for f in frames
                ]
                for ds_id in self.fem.ids
            }
            self._write_fem(w, view, frames[0].fem, (paths, timing))

        self._animate_bodies(w, view, frames, low, high, timing, discrete_timing)

        # the contacts of each step, shown at their step only
        if draw and total:
            arrows, chains, cones = self._contact_layers(most)
            for k, frame in enumerate(frames):
                if len(frame.contacts) == 0:
                    continue
                # visible from its key time to the next one: three values,
                # not one per frame
                values, keys = [], []
                if k > 0:
                    values.append("hidden")
                    keys.append(0.0)
                values.append("visible")
                keys.append(key_times[k] if k > 0 else 0.0)
                if k + 1 < n:
                    values.append("hidden")
                    keys.append(key_times[k + 1])
                extra = (
                    f' visibility="{values[0]}"><animate attributeName="visibility" '
                    f'values="{";".join(values)}" '
                    f'keyTimes="{";".join(f"{t:.6f}" for t in keys)}" '
                    f'dur="{duration:.6g}s" calcMode="discrete" repeatCount="indefinite"/'
                )
                follow = None
                if not discrete and k + 1 < n:
                    follow = self._contacts_follow(
                        view,
                        frame,
                        frames[k + 1],
                        key_times[k],
                        key_times[k + 1],
                        duration,
                    )
                self._write_contacts(
                    w,
                    view,
                    frame,
                    arrows,
                    chains,
                    cones,
                    f"contacts{k}",
                    extra,
                    follow,
                )
        self._write_legend(w, view, low, high)

        # a time bar
        y = view.height - 4
        w(
            f'<g id="clock" font-size="11" fill="{INK}">'
            f'<rect x="8" y="{y - 3:.0f}" width="{view.width - 16:.0f}" height="3" '
            f'fill="{STATIC_FILL}"/>'
            f'<rect x="8" y="{y - 3:.0f}" width="0" height="3" fill="{MUTED}">'
            f'<animate attributeName="width" values="0;{view.width - 16:.0f}" '
            f'dur="{duration:.6g}s" repeatCount="indefinite"/></rect></g>\n'
        )
        w("</svg>\n")

    @staticmethod
    def _contacts_follow(view, frame, after, start, end, duration):
        """follow(ds1) of _write_contacts for the contacts of frame, shown from
        key time start to end: they move with their body ds1 as the body is
        interpolated to its pose in the frame after (a translation of its
        centre, a rotation about it), not left behind it."""

        def pose(f, body):
            k = int(np.searchsorted(f.ids, body))
            if k == len(f.ids) or f.ids[k] != body:
                return None
            c = f.centres[k]
            return (
                (c[0] - view.x0) * view.scale,
                (view.y1 - c[1]) * view.scale,
                f.angles[k],
            )

        def follow(body):
            p0, p1 = pose(frame, body), pose(after, body)
            if p0 is None or p1 is None:
                return "", ""
            (x0, y0, a0), (x1, y1, a1) = p0, p1
            # the rotation of the body between the two frames, as its
            # unwrapped rotation, in svg degrees (y down)
            turn = -np.degrees(np.angle(np.exp(1j * (a1 - a0))))
            keys = [0.0] + ([start] if start > 0 else []) + [end]
            moved = 1 + (start > 0)
            if end < 1:
                keys.append(1.0)
            timing = (
                f' keyTimes="{";".join(f"{t:.6f}" for t in keys)}" '
                f'dur="{duration:.6g}s" repeatCount="indefinite"'
            )
            shift = ";".join(
                "0.00,0.00" if j < moved else f"{x1 - x0:.2f},{y1 - y0:.2f}"
                for j in range(len(keys))
            )
            angle = ";".join(
                f"{0.0 if j < moved else turn:.3f} {x0:.2f} {y0:.2f}"
                for j in range(len(keys))
            )
            opening = (
                f'<g><animateTransform attributeName="transform" type="translate" '
                f'values="{shift}"{timing}/>'
                f'<g><animateTransform attributeName="transform" type="rotate" '
                f'values="{angle}"{timing}/>\n'
            )
            return opening, "</g></g>\n"

        return follow

    def _animate_bodies(self, w, view, frames, low, high, timing, discrete_timing):
        """One group per body: a translation and a rotation animated over the
        frames, the fill too if coloured; hidden where the body is absent."""
        n = len(frames)
        ids = np.unique(np.concatenate([f.ids for f in frames]))
        k_of = {int(i): k for k, i in enumerate(ids)}
        x = np.full((len(ids), n), np.nan)
        y = np.full((len(ids), n), np.nan)
        a = np.zeros((len(ids), n))
        fills = np.full((len(ids), n), BODY, dtype=object)
        present = np.zeros((len(ids), n), dtype=bool)
        block = np.zeros(len(ids), dtype=np.int64)
        scale = np.ones(len(ids))
        for j, frame in enumerate(frames):
            rows = np.array([k_of[int(i)] for i in frame.ids], dtype=np.int64)
            if len(rows) == 0:
                continue
            x[rows, j], y[rows, j] = frame.centres[:, 0], frame.centres[:, 1]
            a[rows, j] = frame.angles
            fills[rows, j] = self._fills(frame.values, low, high)
            present[rows, j] = True
            block[rows] = frame.data["shape_index"]
            scale[rows] = frame.data["glyph_scale"]
        w(f'<g id="bodies" stroke="{INK}" stroke-width="0.8" fill-opacity="0.85">\n')
        for k, i in enumerate(ids):
            mask = present[k]
            first = int(np.argmax(mask))
            # the poses where absent: the nearest earlier one, else the first
            for j in range(n):
                if not mask[j]:
                    src = j - 1 if j > 0 else first
                    x[k, j], y[k, j], a[k, j] = x[k, src], y[k, src], a[k, src]
                    fills[k, j] = fills[k, src]
            if not mask[0]:
                x[k, 0], y[k, 0], a[k, 0] = x[k, first], y[k, first], a[k, first]
                fills[k, 0] = fills[k, first]
            # no jump of 2 pi in the interpolated rotation
            degrees = -np.degrees(np.unwrap(a[k]))
            positions = [view.xy((x[k, j], y[k, j])) for j in range(n)]
            w(f'<g id="b{i}" transform="translate({positions[0]})">')
            if not mask.all():
                w(
                    '<animate attributeName="visibility" values="'
                    + ";".join("visible" if m else "hidden" for m in mask)
                    + f'"{discrete_timing}/>'
                )
            w(
                f'<animateTransform attributeName="transform" type="translate" '
                f'values="{";".join(positions)}"{timing}/>'
            )
            fill = fills[k]
            w(f'<g transform="rotate({degrees[0]:.3f})" fill="{fill[0]}">')
            w(
                f'<animateTransform attributeName="transform" type="rotate" '
                f'values="{";".join(f"{d:.3f}" for d in degrees)}"{timing}/>'
            )
            if len(set(fill)) > 1:
                w(
                    f'<animate attributeName="fill" values="{";".join(fill)}"'
                    f"{discrete_timing}/>"
                )
            self._shape(w, view, block[k], scale[k])
            w("</g></g>\n")
        w("</g>\n")


def write_coulomb(out, step, scheme=None, length_unit=1.0):
    """The Coulomb plot of a step: r_t / (mu r_n) of its loaded contacts
    against r_n on a log scale (the impulses of a pile span decades), the
    cone edges r_t = +- mu r_n as the lines +- 1, the fraction of the
    contacts with r_n > 0 at the limit (|r_t| >= (1 - SATURATION) mu r_n).
    Impulses, N s."""
    contacts = step.contacts
    local = contacts[:, CF["reaction_local"]]
    rn, rt, mu = local[:, 0], local[:, 1], contacts[:, CF["mu"]]
    loaded = rn > 0
    at_limit = loaded & (np.abs(rt) >= (1 - SATURATION) * mu * rn)
    mus = np.unique(np.round(mu, 6)) if len(mu) else np.zeros(0)
    # the frictionless contacts have no ratio: counted, not drawn
    drawn = loaded & (mu > 0)
    ratio = np.zeros(len(rn))
    ratio[drawn] = rt[drawn] / (mu[drawn] * rn[drawn])
    width, height, margin = 480.0, 360.0, 48.0
    if drawn.any():
        lo = int(np.floor(np.log10(rn[drawn].min())))
        hi = int(np.ceil(np.log10(rn[drawn].max())))
    else:
        lo, hi = 0, 1
    hi = max(hi, lo + 1)
    ymax = 1.1 * max(1.0, float(np.abs(ratio).max()) if len(ratio) else 1.0)

    def px(n, t):
        u = (np.log10(n) - lo) / (hi - lo)
        return (
            f"{margin + u * (width - 2 * margin):.2f},"
            f"{height / 2 - t / ymax * (height / 2 - margin):.2f}"
        )

    w = out.write
    w('<?xml version="1.0" encoding="UTF-8"?>\n')
    w(
        f'<svg xmlns="http://www.w3.org/2000/svg" width="{width:.0f}" '
        f'height="{height:.0f}" viewBox="0 0 {width:.0f} {height:.0f}" '
        f'font-family="sans-serif" font-size="11">\n'
    )
    w(f'<rect width="100%" height="100%" fill="{BACKGROUND}"/>\n')
    # the axes, a tick per decade (at most 8 labels)
    left, right = 10.0**lo, 10.0**hi
    w(
        f'<path d="M{px(left, -ymax)}L{px(left, ymax)}M{px(left, 0)}L{px(right, 0)}" '
        f'stroke="{MUTED}" fill="none"/>\n'
    )
    every = -(-(hi - lo) // 8)
    for e in range(lo, hi + 1, every):
        x, y = px(10.0**e, 0).split(",")
        w(
            f'<path d="M{x},{float(y) - 3:.2f}v6" stroke="{MUTED}"/>'
            f'<text x="{x}" y="{float(y) + 14:.2f}" text-anchor="middle" '
            f'fill="{MUTED}">1e{e}</text>\n'
        )
    w(
        f'<text x="{width - margin:.0f}" y="{height - 28:.0f}" text-anchor="end" fill="{INK}">r_n ({_units(length_unit)[1]}), log scale</text>\n'
    )
    for edge in (1, -1):
        x, y = px(left, edge).split(",")
        w(
            f'<text x="{float(x) - 4:.2f}" y="{float(y) + 4:.2f}" text-anchor="end" '
            f'fill="{MUTED}">{edge:+d}</text>\n'
        )
    w(
        f'<text x="{margin + 4:.0f}" y="{margin - 6:.0f}" fill="{INK}">r_t / (mu r_n)</text>\n'
    )
    # the cone edges
    for edge in (1.0, -1.0):
        w(
            f'<path d="M{px(left, edge)}L{px(right, edge)}" '
            f'stroke="{CONE}" stroke-dasharray="4,3" fill="none"/>\n'
        )
    # the contacts, those at the limit in red
    for selected, colour in ((drawn & ~at_limit, INK), (drawn & at_limit, FORCE)):
        if selected.any():
            d = "".join(f"M{px(a, b)}h0" for a, b in zip(rn[selected], ratio[selected]))
            w(
                f'<path d="{d}" stroke="{colour}" stroke-width="4" '
                f'stroke-linecap="round" fill="none"/>\n'
            )
    frictionless = int((loaded & ~drawn).sum())
    if frictionless:
        w(
            f'<text x="{width - margin:.0f}" y="{margin - 6:.0f}" text-anchor="end" '
            f'fill="{INK}">{frictionless} frictionless, not drawn</text>\n'
        )
    n_loaded, n_limit = int(loaded.sum()), int(at_limit.sum())
    fraction = f"{100 * n_limit / n_loaded:.1f}%" if n_loaded else "-"
    title = f"t = {step.time:.6g} s" + (f" · {scheme}" if scheme else "")
    title += f" · mu = {', '.join(f'{m:.3g}' for m in mus)}" if len(mus) else ""
    w(f'<text x="8" y="16" font-size="12" fill="{INK}">{_escape(title)}</text>\n')
    w(
        f'<text x="8" y="{height - 10:.0f}" fill="{INK}">at the limit '
        f"(|r_t| ≥ (1 - {SATURATION:g}) mu r_n): {n_limit} of {n_loaded} "
        f"loaded contacts, {fraction}</text>\n"
    )
    w("</svg>\n")
    return n_limit, n_loaded


def coulomb_filename(io_filename, time, directory=None):
    """<stem>_coulomb_<time>.svg, as svg_filename."""
    name = svg_filename(io_filename, time, directory)
    return name[: -len(f"_svg_{time:013.6f}.svg")] + f"_coulomb_{time:013.6f}.svg"


def animation_filename(io_filename, directory=None, kind="animation"):
    """<stem>_svg_animation.svg (or _svg_strobe.svg) beside the HDF5 file."""
    stem = os.path.splitext(os.path.basename(io_filename))[0]
    directory = directory if directory is not None else os.path.dirname(io_filename)
    return os.path.join(directory, f"{stem}_svg_{kind}.svg")


# -- during a run ---------------------------------------------------------------


class RunPictures:
    """The pictures of the svg_output and svg_options run options.

    requests: times (float) and simulation step numbers (int, the time k h).
    A request is drawn at the step k with t_k <= request < t_k + h, an output
    step or not. The votes for the sides are those of the output steps and of
    the drawn steps, as convert counts those of the file.
    options: those of Picture, and coulomb (the Coulomb plot of each
    picture), animate and strobe (of the output steps, written at the end;
    the steps are held in memory), uniform, discrete, speed (of the
    animation), directory.
    """

    def __init__(self, io_filename, shapes, requests, h, fem=None, **options):
        self.io_filename = io_filename
        self.directory = options.pop("directory", None)
        self.coulomb = options.pop("coulomb", False)
        self.animate = options.pop("animate", False)
        self.strobe = options.pop("strobe", False)
        self.animation_options = {
            key: options.pop(key)
            for key in ("uniform", "discrete", "speed")
            if key in options
        }
        self.h = h
        times = [
            k * h if isinstance(k, (int, np.integer)) else float(k) for k in requests
        ]
        self.pending = sorted(times)
        self.picture = Picture(shapes, fem, **options)
        self._voted = None
        self.steps = []
        self.names = []

    def preload(self, data):
        """A run continued from its file: the votes of the output steps of
        the file, as convert counts them."""
        if "step_index" not in data:
            return
        steps = FileSteps(data)
        for k in range(len(steps.times)):
            step = steps.step(k)
            if len(step.bodies):
                self.picture.first_bodies(step.bodies)
            self.picture.vote(step)
            self._voted = step.time

    def due(self, time):
        """The requests of the step at time."""
        eps = 1e-9 * self.h
        return [t for t in self.pending if time - eps <= t < time + self.h - eps]

    def output_step(self, step):
        """Count the votes of an output step; the first one gives the bodies of
        the segments without contact. Keep it for the animation."""
        self.picture.first_bodies(step.bodies)
        if self._voted != step.time:
            self.picture.vote(step)
            self._voted = step.time
            if self.animate or self.strobe:
                self.steps.append(step)

    def draw(self, step):
        """The pictures due at the step (rows as the tables have them)."""
        for requested in self.due(step.time):
            self.picture.first_bodies(step.bodies)
            if self._voted != step.time:
                self.picture.vote(step)
                self._voted = step.time
            name = svg_filename(self.io_filename, step.time, self.directory)
            with open(name, "w") as out:
                self.picture.write(out, step, requested)
            self.names.append(name)
            if self.coulomb:
                name = coulomb_filename(self.io_filename, step.time, self.directory)
                with open(name, "w") as out:
                    write_coulomb(
                        out, step, self.picture.scheme, self.picture.length_unit
                    )
                self.names.append(name)
            if abs(requested - step.time) > 1e-9 * self.h:
                warn(
                    f"time {requested:.6g} requested, the step at {step.time:.6g} drawn"
                )
            self.pending.remove(requested)

    def finish(self):
        for requested in self.pending:
            warn(f"time {requested:.6g} requested, not reached by the run: not drawn")
        for kind, wanted in (("animation", self.animate), ("strobe", self.strobe)):
            if wanted and self.steps:
                name = animation_filename(self.io_filename, self.directory, kind)
                with open(name, "w") as out:
                    if kind == "animation":
                        self.picture.write_animation(
                            out, self.steps, **self.animation_options
                        )
                    else:
                        self.picture.write_strobe(out, self.steps)
                self.names.append(name)


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
        # the FEM tables: ds id -> (times, table)
        self.fem = {}
        for name in data:
            if name.startswith("fem_displacements_") and data[name].shape[0]:
                table = data[name]
                self.fem[int(name[len("fem_displacements_") :])] = (
                    np.asarray(table[:, 0], dtype=np.float64),
                    table,
                )

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
        # a time equal to an output step up to the rounding is that step
        eps = SAME_TIME * max(1.0, abs(time))
        k = int(np.searchsorted(self.times, time + eps, side="right")) - 1
        if k < 0:
            raise ValueError(
                f"time {time} before the first output step, {self.times[0]}"
            )
        if not _same_time(time, self.times[k]):
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

    def fem_rows(self, time):
        """{ds id: displacements} of the FEM tables at time (none before their
        first row: undeformed)."""
        rows = {}
        for ds_id, (times, table) in self.fem.items():
            k = np.flatnonzero(times == time)
            if len(k):
                rows[ds_id] = np.asarray(table[int(k[0]), 2:], dtype=np.float64)
        return rows

    def step(self, k):
        """The Step of output step k."""
        time = float(self.times[k])
        return Step(
            time,
            self.rows("dynamic", k, 8),
            self.rows("velocities", k, 7),
            self.rows("static", k, 8),
            self.rows("cf", k, self.contact_columns),
            self.fem_rows(time),
        )


def convert(
    io_filename,
    times=(),
    frames=(),
    directory=None,
    output=None,
    coulomb=False,
    animate=False,
    strobe=False,
    stride=1,
    uniform=False,
    discrete=False,
    speed=1.0,
    **options,
):
    """Write the pictures of the output steps at times (at or before each)
    and frames (output step indices) of a siconos mechanics HDF5 file, the
    last step if none (and no animation).

    coulomb: the Coulomb plot of each picture too.
    animate, strobe: an animation, a strobe picture of the output steps
    (every stride-th one, the last included); the pictures of times and
    frames then use the view of the animation, so that a picture is the
    frame of the animation.
    options: those of Picture. Returns the file names.
    """
    import h5py

    if directory is not None:
        os.makedirs(directory, exist_ok=True)
    names = []
    with h5py.File(io_filename, "r") as f:
        dimension = int(f.attrs.get("dimension", 3))
        if dimension != 2:
            raise ValueError(f"{io_filename}: SVG pictures are for 2D files")
        data = f["data"]
        steps = FileSteps(data)
        run_options = _file_run_options(data)
        options.setdefault("scheme", run_options.get("osi"))
        options.setdefault("length_unit", run_options.get("gravity_scale"))
        picture = Picture(
            vtkhdf.ShapeTable(data, dimension), vtkhdf.FemMeshes(data), **options
        )
        n = len(steps.times)
        selected = sorted(set(range(0, n, max(stride, 1))) | {n - 1})
        requests = [(steps.resolve(time=t), t) for t in times]
        if animate or strobe:
            # a frame of the animation
            for k in frames:
                if not -len(selected) <= k < len(selected):
                    raise ValueError(f"frame {k} not in 0..{len(selected) - 1}")
                requests.append((selected[k], None))
        else:
            requests += [(steps.resolve(frame=k), None) for k in frames]
        if not requests and not (animate or strobe):
            requests = [(steps.resolve(), None)]
        if output is not None and len(requests) + bool(animate) + bool(strobe) > 1:
            raise ValueError("one output file name for several pictures")
        picture.first_bodies(steps.step(0).bodies)
        bounds = None
        if animate or strobe:
            # all the votes: the sides of an animation do not change
            for k in range(n):
                picture.vote(steps.step(k))
            animation_steps = [steps.step(k) for k in selected]
            bounds = picture.bounds([_Frame(picture, s) for s in animation_steps])
            for kind, wanted in (("animation", animate), ("strobe", strobe)):
                if not wanted:
                    continue
                name = output or animation_filename(io_filename, directory, kind)
                with open(name, "w") as out:
                    if kind == "animation":
                        picture.write_animation(
                            out, animation_steps, uniform, discrete, speed
                        )
                    else:
                        picture.write_strobe(out, animation_steps, bounds)
                names.append(name)
            voted = n - 1
        else:
            voted = -1
        # the votes of the steps up to each drawn one, in order
        for k, requested in sorted(requests, key=lambda r: r[0]):
            for j in range(voted + 1, k + 1):
                picture.vote(steps.step(j))
            voted = max(voted, k)
            step = steps.step(k)
            name = output or svg_filename(io_filename, step.time, directory)
            with open(name, "w") as out:
                picture.write(out, step, requested, bounds)
            names.append(name)
            if coulomb:
                name = coulomb_filename(io_filename, step.time, directory)
                with open(name, "w") as out:
                    write_coulomb(out, step, picture.scheme, picture.length_unit)
                names.append(name)
    return names
