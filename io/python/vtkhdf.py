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
"""Transient VTKHDF output of the bodies of a mechanics simulation.

The bodies are written as one transient PolyData, one point per body,
read natively by ParaView (>= 5.12) and VTK (>= 9.3). It is meant for
large granular simulations. Alongside <name>.vtkhdf:

- <name>_shapes.vtkhdf, the glyph table: one block per distinct shape of
  the bodies in the body frame, picked for each body by its shape_index
  point array, scaled by glyph_scale (the radius of the disks and
  spheres);
- <name>_static.vtkhdf, the static objects alive, their geometry merged
  in one PolyData, a step written only when it changes;
- <name>_contacts.vtkhdf, the contacts of each step, one point per
  contact (the middle of its points A and B) with its normal and reaction;
- <name>_fem.vtkhdf, the finite element meshes deformed at each step,
  with their displacement (point data), strain, stress and von Mises
  stress (cell data);
- <name>_paraview.py, a ParaView state drawing every body and static
  object with its shape (3D Glyphs oriented by the "orientation"
  quaternion w, x, y, z): paraview --state=<name>_paraview.py

For the large cases, BODIES_AS_DISKS = True at the top of the state (the
disks option of the writers) draws every body as a disk scaled by the
"disk_radius" array, with the "Point Gaussian" representation: the
3D Glyphs of ParaView keep a copy of the glyph table at every step shown.

Layout: every step appends its rows to the same datasets, and the
Steps group records where each step starts. The vertex topology (one
vertex per point) is written once, as the range 0..max(N), and every
step reads a prefix of it.

Only numpy and h5py are needed, not VTK.
"""

import os
from collections import namedtuple

import h5py
import numpy as np

VTKHDF_VERSION = (2, 1)

# chunk rows of the datasets with one row per step
_STEP_CHUNK_ROWS = 256

# the topology groups of a PolyData, in the order of the columns of
# Steps/CellOffsets and Steps/ConnectivityIdOffsets
_TOPOLOGIES = ("Vertices", "Lines", "Polygons", "Strips")


def _extendable(group, name, shape, dtype, chunk_rows, compression):
    """A dataset of shape (0, *shape), extendable along the first axis."""
    return group.create_dataset(
        name,
        (0, *shape),
        maxshape=(None, *shape),
        chunks=(chunk_rows, *shape),
        dtype=dtype,
        compression=compression,
        compression_opts=4 if compression == "gzip" else None,
    )


def _append(dataset, values):
    n = dataset.shape[0]
    dataset.resize(n + len(values), axis=0)
    dataset[n:] = values
    return n


class BodiesWriter:
    """Append the state of the bodies, step by step, to a VTKHDF file.

    Each call to write_step adds one time step. The number of bodies
    may change between steps (births and deaths).

    Parameters
    ----------
    filename : str
        the .vtkhdf file, overwritten unless append.
    dtype : numpy dtype, optional
        storage type of the coordinates and of the point arrays,
        float32 by default, like the siconos HDF5 tables.
    compression : bool, optional
        gzip the datasets, default False.
    chunk_rows : int, optional
        number of rows of an HDF5 chunk of the point data. The per step
        datasets (Steps, NumberOfPoints, ...) use small chunks.
    append : bool, optional
        continue the steps of an existing file (a restart), default
        False. A missing file is created.
    """

    def __init__(
        self,
        filename,
        dtype=np.float32,
        compression=False,
        chunk_rows=8192,
        append=False,
    ):
        self._dtype = dtype
        self._compression = "gzip" if compression else None
        self._chunk_rows = chunk_rows
        self._number_of_steps = 0
        self._number_of_points = 0
        self._vertices_size = 0
        self._point_arrays = {}

        if append and os.path.exists(filename):
            self._file = h5py.File(filename, "a")
            self._reopen()
        else:
            self._file = h5py.File(filename, "w")
            self._create()

    def _reopen(self):
        """Take up the state of the steps already in the file."""
        root = self._root = self._file["VTKHDF"]
        self._steps = root["Steps"]
        self._points = root["Points"]
        self._number_of_points_ds = root["NumberOfPoints"]
        self._number_of_steps = int(self._steps.attrs["NSteps"])
        self._number_of_points = self._points.shape[0]
        self._vertices_size = root["Vertices"]["Connectivity"].shape[0]
        self._point_arrays = dict(root["PointData"].items())

    def _create(self):
        dtype = self._dtype
        root = self._file.create_group("VTKHDF")
        root.attrs["Version"] = np.array(VTKHDF_VERSION, dtype=np.int64)
        # the reader wants a fixed length ASCII string, not a variable one
        root.attrs["Type"] = np.bytes_("PolyData")
        self._root = root

        self._points = self._extendable(root, "Points", (3,), dtype)
        self._number_of_points_ds = self._per_step(root, "NumberOfPoints", (), np.int64)

        for name in _TOPOLOGIES:
            topology = root.create_group(name)
            for count in ("NumberOfCells", "NumberOfConnectivityIds"):
                self._per_step(topology, count, (), np.int64)
            if name == "Vertices":
                self._extendable(topology, "Offsets", (), np.int64)
                self._extendable(topology, "Connectivity", (), np.int64)
            else:
                # no cell, but a step still reads one offset
                self._per_step(topology, "Offsets", (), np.int64)
                self._per_step(topology, "Connectivity", (), np.int64)
                _append(topology["Offsets"], [0])

        root.create_group("PointData")
        root.create_group("CellData")
        root.create_group("FieldData")

        steps = root.create_group("Steps")
        steps.attrs["NSteps"] = 0
        self._steps = steps
        self._per_step(steps, "Values", (), np.float64)
        for name in ("PartOffsets", "NumberOfParts", "PointOffsets"):
            self._per_step(steps, name, (), np.int64)
        for name in ("CellOffsets", "ConnectivityIdOffsets"):
            self._per_step(steps, name, (len(_TOPOLOGIES),), np.int64)
        steps.create_group("PointDataOffsets")
        steps.create_group("CellDataOffsets")
        steps.create_group("FieldDataOffsets")

    def _extendable(self, group, name, shape, dtype):
        return _extendable(
            group, name, shape, dtype, self._chunk_rows, self._compression
        )

    def _per_step(self, group, name, shape, dtype):
        return _extendable(group, name, shape, dtype, _STEP_CHUNK_ROWS, None)

    def _grow_vertices(self, n):
        """Make the vertex topology 0..n-1 available as a prefix."""
        if n <= self._vertices_size:
            return
        vertices = self._root["Vertices"]
        first = self._vertices_size
        if first == 0:
            _append(vertices["Offsets"], [0])
        _append(vertices["Offsets"], np.arange(first + 1, n + 1, dtype=np.int64))
        _append(vertices["Connectivity"], np.arange(first, n, dtype=np.int64))
        self._vertices_size = n

    def _point_array(self, name, width, integer):
        """The dataset of a point array, created on its first use."""
        dataset = self._point_arrays.get(name)
        if dataset is None:
            if self._number_of_steps > 0:
                raise ValueError(
                    f"point array {name} must be given from the first step"
                )
            shape = () if width == 1 else (width,)
            dtype = np.int64 if integer else self._dtype
            dataset = self._extendable(self._root["PointData"], name, shape, dtype)
            self._per_step(self._steps["PointDataOffsets"], name, (), np.int64)
            self._point_arrays[name] = dataset
        return dataset

    def write_step(self, time, points, point_data):
        """Append one time step.

        Parameters
        ----------
        time : float
        points : array (N, 2) or (N, 3)
            the body positions, z = 0 in 2D.
        point_data : dict
            name -> array of N rows, the same names at every step.
        """
        points = np.asarray(points)
        n = points.shape[0]
        if points.shape[1] == 2:
            points = np.column_stack((points, np.zeros(n)))

        if self._number_of_steps > 0 and set(point_data) != set(self._point_arrays):
            raise ValueError(
                f"point arrays {sorted(point_data)} differ from the first step "
                f"{sorted(self._point_arrays)}"
            )

        steps = self._steps
        _append(steps["Values"], [time])
        _append(steps["PartOffsets"], [self._number_of_steps])
        _append(steps["NumberOfParts"], [1])
        _append(steps["PointOffsets"], [self._number_of_points])
        # every step reads its vertices from the start of the shared range
        _append(steps["CellOffsets"], np.zeros((1, len(_TOPOLOGIES)), dtype=np.int64))
        _append(
            steps["ConnectivityIdOffsets"],
            np.zeros((1, len(_TOPOLOGIES)), dtype=np.int64),
        )

        _append(self._number_of_points_ds, [n])
        _append(self._points, points.astype(self._dtype, copy=False))

        self._grow_vertices(n)
        for name in _TOPOLOGIES:
            topology = self._root[name]
            cells = n if name == "Vertices" else 0
            _append(topology["NumberOfCells"], [cells])
            _append(topology["NumberOfConnectivityIds"], [cells])

        for name, values in point_data.items():
            values = np.asarray(values)
            if values.shape[0] != n:
                raise ValueError(
                    f"point array {name}: {values.shape[0]} rows, {n} points"
                )
            width = 1 if values.ndim == 1 else values.shape[1]
            integer = name == "id" or np.issubdtype(values.dtype, np.integer)
            dataset = self._point_array(name, width, integer)
            offset = _append(dataset, values.astype(dataset.dtype, copy=False))
            _append(steps["PointDataOffsets"][name], [offset])

        self._number_of_points += n
        self._number_of_steps += 1
        steps.attrs["NSteps"] = self._number_of_steps

    def last_step(self):
        """(time, points, point_data) of the last step, None if no step."""
        k = self._number_of_steps - 1
        if k < 0:
            return None
        steps = self._steps
        n = int(self._number_of_points_ds[k])
        first = int(steps["PointOffsets"][k])
        point_data = {}
        for name, dataset in self._point_arrays.items():
            offset = int(steps["PointDataOffsets"][name][k])
            point_data[name] = dataset[offset : offset + n]
        return float(steps["Values"][k]), self._points[first : first + n], point_data

    def flush(self):
        self._file.flush()

    def close(self):
        if self._file:
            self._file.close()
            self._file = None

    def __enter__(self):
        return self

    def __exit__(self, type_, value, traceback):
        self.close()


# --- shapes ----------------------------------------------------------------
#
# Each body is drawn with the "3D Glyphs" representation and a glyph
# table: the shapes file is a MultiBlockDataSet, one block per distinct
# shape in the body frame, and the shape_index point array picks the block
# of each body. Disks and spheres share a unit block scaled by the
# glyph_scale point array (their radius); the other shapes are at their
# size, glyph_scale = 1.

# number of sides of a disk, of meridians of a sphere or a cylinder
_RESOLUTION = 24

# drawn length of the infinite Line primitive
_LINE_LENGTH = 1000.0


def _unit_circle(n=_RESOLUTION):
    a = 2 * np.pi * np.arange(n) / n
    return np.column_stack((np.cos(a), np.sin(a)))


def _turn(o, a, b):
    """> 0 for a counterclockwise turn o -> a -> b."""
    return (a[0] - o[0]) * (b[1] - o[1]) - (a[1] - o[1]) * (b[0] - o[0])


def _hull_2d(points):
    """Convex hull of 2D points, counterclockwise (monotone chain)."""
    pts = sorted(map(tuple, np.asarray(points, dtype=np.float64)[:, :2]))
    if len(pts) < 3:
        return np.array(pts)

    def half(seq):
        h = []
        for p in seq:
            while len(h) >= 2 and _turn(h[-2], h[-1], p) <= 0:
                h.pop()
            h.append(p)
        return h[:-1]

    return np.array(half(pts) + half(reversed(pts)))


class _Geometry:
    """Points (N, 3) and polygonal cells, the cells as lists of point ids."""

    def __init__(self, points=None, cells=()):
        self.points = np.zeros((0, 3)) if points is None else _xyz(points)
        self.cells = [np.asarray(c, dtype=np.int64) for c in cells]

    def merge(self, other, translation=(0, 0, 0), orientation=(1, 0, 0, 0)):
        n = len(self.points)
        moved = _rotate(orientation, other.points) + np.asarray(translation)[:3]
        self.points = np.vstack((self.points, moved))
        self.cells += [c + n for c in other.cells]
        return self


def _xyz(points):
    points = np.asarray(points, dtype=np.float64)
    if points.ndim == 2 and points.shape[1] == 2:
        points = np.column_stack((points, np.zeros(len(points))))
    return points


def _rotate(q, points):
    """Rotate points (N, 3) by the unit quaternion q = (w, x, y, z)."""
    w, v = q[0], np.asarray(q[1:4], dtype=np.float64)
    return points + 2 * np.cross(v, np.cross(v, points) + w * points)


def _polygon_area(points):
    """Area of a planar polygon of the xy plane (shoelace)."""
    x, y = points[:, 0], points[:, 1]
    return 0.5 * abs(np.dot(x, np.roll(y, -1)) - np.dot(y, np.roll(x, -1)))


def _polygon(points):
    return _Geometry(points, [np.arange(len(points))])


def _sphere(radius=1.0, n=_RESOLUTION):
    """A UV sphere, triangles at the poles, quads elsewhere."""
    rings = n // 2
    theta = np.pi * np.arange(1, rings) / rings
    circle = _unit_circle(n)
    points = [(0, 0, radius)]
    for t in theta:
        points += [
            (radius * np.sin(t) * c, radius * np.sin(t) * s, radius * np.cos(t))
            for c, s in circle
        ]
    points.append((0, 0, -radius))
    south = len(points) - 1
    cells = [(0, 1 + i, 1 + (i + 1) % n) for i in range(n)]
    for r in range(rings - 2):
        a, b = 1 + r * n, 1 + (r + 1) * n
        cells += [(a + i, b + i, b + (i + 1) % n, a + (i + 1) % n) for i in range(n)]
    last = 1 + (rings - 2) * n
    cells += [(last + i, south, last + (i + 1) % n) for i in range(n)]
    return _Geometry(np.array(points), cells)


def _box(lx, ly, lz):
    corners = np.array([(x, y, z) for z in (-1, 1) for y in (-1, 1) for x in (-1, 1)])
    faces = [
        (0, 2, 3, 1),
        (4, 5, 7, 6),
        (0, 1, 5, 4),
        (2, 6, 7, 3),
        (0, 4, 6, 2),
        (1, 3, 7, 5),
    ]
    return _Geometry(corners * np.array([lx, ly, lz]) / 2, faces)


def _cylinder(radius, height, top_radius=None, n=_RESOLUTION):
    """Along y, centered, like vtkCylinderSource; a cone if top_radius = 0."""
    top_radius = radius if top_radius is None else top_radius
    circle = _unit_circle(n)
    bottom = np.column_stack(
        (radius * circle[:, 0], np.full(n, -height / 2), radius * circle[:, 1])
    )
    if top_radius == 0:
        points = np.vstack((bottom, [(0, height / 2, 0)]))
        cells = [(i, (i + 1) % n, n) for i in range(n)] + [np.arange(n)[::-1]]
        return _Geometry(points, cells)
    top = np.column_stack(
        (top_radius * circle[:, 0], np.full(n, height / 2), top_radius * circle[:, 1])
    )
    cells = [(i, (i + 1) % n, n + (i + 1) % n, n + i) for i in range(n)]
    cells += [np.arange(n)[::-1], n + np.arange(n)]
    return _Geometry(np.vstack((bottom, top)), cells)


def _convex_3d(vertices):
    from scipy.spatial import ConvexHull

    hull = ConvexHull(vertices)
    return _Geometry(vertices, hull.simplices)


def _shape_geometry(shape, dimension):
    """(geometry, unit): the geometry of a shape of the ref group, unit for
    the disks and spheres drawn as a unit block scaled by their radius.
    None if the shape is not drawn."""
    kind = shape.attrs.get("type")
    if isinstance(kind, bytes):
        kind = kind.decode()
    if kind == "convex":
        vertices = np.asarray(shape[:], dtype=np.float64)
        if vertices.shape[1] == 2:
            return _polygon(_hull_2d(vertices)), None
        return _convex_3d(vertices), None
    if kind != "primitive":
        return None
    primitive = shape.attrs["primitive"]
    if isinstance(primitive, bytes):
        primitive = primitive.decode()
    params = np.asarray(shape[:], dtype=np.float64).ravel()
    if primitive in ("Disk", "Circle"):
        return None, ("disk", params[0])
    if primitive == "Sphere":
        return None, ("sphere", params[0])
    if primitive == "Segment":
        x1, y1, x2, y2 = params[:4]
        return _Geometry(np.array([(x1, y1), (x2, y2)]), [(0, 1)]), None
    if primitive == "Line":
        # a x + b y + c = 0, drawn over _LINE_LENGTH around its nearest
        # point to the origin
        a, b, c = params[:3]
        n2 = a * a + b * b
        p0 = np.array((-a * c / n2, -b * c / n2))
        d = np.array((-b, a)) / np.sqrt(n2) * _LINE_LENGTH / 2
        return _Geometry(np.array([p0 - d, p0 + d]), [(0, 1)]), None
    if primitive == "Box2d":
        lx, ly = params[:2] / 2
        return _polygon([(-lx, -ly), (lx, -ly), (lx, ly), (-lx, ly)]), None
    if primitive == "Box":
        return _box(*params[:3]), None
    if primitive == "Cylinder":
        return _cylinder(params[0], params[1]), None
    if primitive == "Cone":
        return _cylinder(params[0], params[1], top_radius=0), None
    if primitive == "Capsule":
        # as vview draws it
        r, h = params[:2]
        g = _cylinder(r, h / 2)
        g.merge(_sphere(r), (0, h / 4, 0))
        return g.merge(_sphere(r), (0, -h / 4, 0)), None
    return None


def _decode(names):
    return [n.decode("utf-8") if isinstance(n, bytes) else str(n) for n in names]


def _cross():
    """A small cross, for the bodies of unknown or unsupported shapes."""
    r = 0.5
    arms = [(-r, 0, 0), (r, 0, 0), (0, -r, 0), (0, r, 0), (0, 0, -r), (0, 0, r)]
    return _Geometry(np.array(arms), [(0, 1), (2, 3), (4, 5)])


def _contactor_offset(c):
    t = np.asarray(c.attrs.get("translation", (0, 0, 0)), dtype=np.float64)
    q = np.asarray(c.attrs.get("orientation", (1, 0, 0, 0)), dtype=np.float64)
    return np.pad(t, (0, 3 - len(t))), q


def _contactors_geometry(contactors, ref, dimension, unsupported):
    """The contactors of an object merged in its frame, at their size; the
    names of the shapes not drawn are added to unsupported."""
    g = _Geometry()
    for c in contactors:
        shape_name = c.attrs["shape_name"]
        part = _shape_geometry(ref[shape_name], dimension)
        if part is None:
            unsupported.add(shape_name)
            continue
        part, unit = part
        if unit is not None:
            r = unit[1]
            part = _polygon(r * _unit_circle()) if unit[0] == "disk" else _sphere(r)
        g.merge(part, *_contactor_offset(c))
    return g


class ShapeTable:
    """The shapes of the bodies of a siconos mechanics HDF5 file.

    Built from the data/input and data/ref groups. The grain k of an
    add_objects aggregate has the id of the aggregate + k.

    blocks: list of (name, _Geometry), the glyph table of the dynamic
    bodies only.
    point_data(ids): the shape_index, glyph_scale and disk_radius
    arrays of bodies.
    static_geometry(rows, time): the static objects (negative ids) alive
    at time, from their time of birth to their time of death, merged in
    one geometry: they are not glyphs (see _PARAVIEW_STATE).
    """

    def __init__(self, data, dimension):
        self.dimension = dimension
        self.blocks = []
        self._block_of_key = {}
        self.unsupported = set()
        self._fallback = None
        # static id -> (geometry in the object frame, time of birth, time
        # of death)
        self._static = {}
        # of the static objects and of the bodies at their initial positions
        self._bounds = None
        ref = data["ref"]
        # id -> (block, scale), filled in dense arrays at the end
        ids, blocks, scales = [], [], []

        def unit_block():
            name = "Disk" if dimension == 2 else "Sphere"
            geometry = _polygon(_unit_circle()) if dimension == 2 else _sphere()
            return self._block(("unit",), name, lambda: geometry)

        def shape_block(name):
            """(block, scale) of a shape alone, at the body origin."""
            g = _shape_geometry(ref[name], dimension)
            if g is None:
                self.unsupported.add(name)
                return self._fallback_block(), 1.0
            geometry, unit = g
            if unit is not None:
                return unit_block(), unit[1]
            return self._block(("shape", name), name, lambda: geometry), 1.0

        for name, obj in data["input"].items():
            kind = obj.attrs.get("type", "")
            if isinstance(kind, bytes):
                kind = kind.decode()
            base = int(obj.attrs["id"])
            if "material" in obj.attrs:
                # a finite element body, drawn from its mesh (FemMeshes)
                continue
            if kind == "dynamic_aggregate":
                count = int(obj.attrs["count"])
                if "shape_names" in obj:
                    names = _decode(obj["shape_names"][:])
                    indices = np.asarray(obj["shape_indices"][:], dtype=np.int64)
                else:
                    names, indices = (
                        [obj.attrs["shape_name"]],
                        np.zeros(count, np.int64),
                    )
                per_shape = [shape_block(n) for n in names]
                block = np.array([b for b, _ in per_shape])[indices]
                scale = np.array([s for _, s in per_shape])[indices]
                unit = block == self._block_of_key.get(("unit",), -1)
                if unit.any() and ("radii" in obj or "radius" in obj.attrs):
                    radii = obj["radii"][:] if "radii" in obj else obj.attrs["radius"]
                    scale = np.where(unit, np.broadcast_to(radii, count), scale)
                ids.append(base + np.arange(count))
                blocks.append(block)
                scales.append(scale)
                self._extend_bounds(
                    obj["translations"][:], self._block_radius(block, scale)
                )
            elif kind == "dynamic" or (kind != "static" and "mass" in obj.attrs):
                contactors = [obj[c] for c in obj if "shape_name" in obj[c].attrs]
                block, scale = self._contactors_block(
                    name, contactors, shape_block, ref
                )
                ids.append([base])
                blocks.append([block])
                scales.append([scale])
                self._extend_bounds(
                    [obj.attrs.get("translation", (0, 0, 0))],
                    self._block_radius([block], [scale]),
                )
            elif kind == "static":
                contactors = [obj[c] for c in obj if "shape_name" in obj[c].attrs]
                geometry = _contactors_geometry(
                    contactors, ref, dimension, self.unsupported
                )
                if len(geometry.points) == 0:
                    geometry = _cross()
                self._static[base] = (
                    geometry,
                    obj.attrs.get("time_of_birth", -np.inf),
                    obj.attrs.get("time_of_death", np.inf),
                )
                orientation = np.asarray(
                    obj.attrs.get("orientation", (1, 0, 0, 0)), dtype=np.float64
                )
                if orientation.size == 1:
                    # a 2D angle
                    half = orientation[0] / 2
                    orientation = np.array((np.cos(half), 0, 0, np.sin(half)))
                given = np.asarray(obj.attrs.get("translation", ()), dtype=np.float64)
                translation = np.zeros(3)
                translation[: min(3, given.size)] = given.ravel()[:3]
                self._extend_bounds(
                    _rotate(orientation, geometry.points) + translation, 0.0
                )

        ids = np.concatenate(ids).astype(np.int64) if ids else np.zeros(0, np.int64)
        size = int(ids.max()) + 1 if len(ids) else 0
        fallback = self._fallback_block()
        self._shape_index = np.full(size, fallback, dtype=np.int64)
        self._glyph_scale = np.ones(size)
        if len(ids):
            self._shape_index[ids] = np.concatenate(blocks)
            self._glyph_scale[ids] = np.concatenate(scales)
        # the largest distance of a block to the body origin
        self._extent = np.array(
            [
                np.linalg.norm(g.points, axis=1).max() if len(g.points) else 0.0
                for _, g in self.blocks
            ]
        )
        # the radius of a block drawn as a disk: in 2D, the disk of the
        # same area as its polygons, else the largest distance
        self._disk_radius = self._extent.copy()
        if dimension == 2:
            for k, (_, g) in enumerate(self.blocks):
                area = sum(_polygon_area(g.points[c]) for c in g.cells if len(c) > 2)
                if area > 0:
                    self._disk_radius[k] = np.sqrt(area / np.pi)
        unit = self._block_of_key.get(("unit",))
        if unit is not None:
            self._disk_radius[unit] = 1.0
        # the median size of the bodies, for the contacts drawing
        self.typical_size = 1.0
        if len(ids):
            sizes = self._extent[self._shape_index[ids]] * self._glyph_scale[ids]
            if np.any(sizes > 0):
                self.typical_size = float(np.median(sizes[sizes > 0]))

    def _block_radius(self, blocks, scales):
        """The largest extent of the blocks from the body origin."""
        extent = 0.0
        for block in np.unique(blocks):
            points = self.blocks[block][1].points
            if len(points):
                largest = np.max(np.asarray(scales)[np.asarray(blocks) == block])
                extent = max(extent, largest * np.linalg.norm(points, axis=1).max())
        return extent

    def _extend_bounds(self, points, margin):
        points = np.asarray(points, dtype=np.float64)
        if points.size == 0:
            return
        points = np.pad(points, ((0, 0), (0, max(0, 3 - points.shape[1]))))[:, :3]
        low, high = points.min(axis=0) - margin, points.max(axis=0) + margin
        if self._bounds is not None:
            low = np.minimum(low, self._bounds[0])
            high = np.maximum(high, self._bounds[1])
        self._bounds = (low, high)

    @property
    def bounds(self):
        """xmin, xmax, ymin, ymax, zmin, zmax of the scene, None if empty."""
        if self._bounds is None:
            return None
        low, high = self._bounds
        return [float(v) for pair in zip(low, high) for v in pair]

    def _block(self, key, name, geometry):
        block = self._block_of_key.get(key)
        if block is None:
            block = len(self.blocks)
            self.blocks.append((name, geometry()))
            self._block_of_key[key] = block
        return block

    def _fallback_block(self):
        """A small cross, for the bodies of unknown or unsupported shapes."""
        if self._fallback is None:
            self._fallback = self._block(("fallback",), "unknown", _cross)
        return self._fallback

    def _contactors_block(self, name, contactors, shape_block, ref):
        """(block, scale) of a body made of its contactors."""

        def identity(c):
            t, q = _contactor_offset(c)
            return np.allclose(t, 0) and len(q) == 4 and np.allclose(q, (1, 0, 0, 0))

        if not contactors:
            return self._fallback_block(), 1.0
        if len(contactors) == 1 and identity(contactors[0]):
            return shape_block(contactors[0].attrs["shape_name"])
        # a composite body: the contactors merged, at their size
        key = tuple(
            (c.attrs["shape_name"], *_contactor_offset(c)[0], *_contactor_offset(c)[1])
            for c in contactors
        )

        def geometry():
            return _contactors_geometry(
                contactors, ref, self.dimension, self.unsupported
            )

        return self._block(("composite", key), name, geometry), 1.0

    def point_data(self, ids):
        """shape_index, glyph_scale and disk_radius (for the bodies drawn
        as disks: in 2D the radius of the disk of the same area as the
        shape, else the largest distance of the shape to the body origin)
        of the bodies ids."""
        ids = np.asarray(ids).astype(np.int64)
        known = (ids >= 0) & (ids < len(self._shape_index))
        index = np.full(len(ids), self._fallback_block(), dtype=np.int64)
        scale = np.ones(len(ids))
        index[known] = self._shape_index[ids[known]]
        scale[known] = self._glyph_scale[ids[known]]
        radius = self._disk_radius[index] * scale
        return {"shape_index": index, "glyph_scale": scale, "disk_radius": radius}

    @property
    def has_static(self):
        return len(self._static) > 0

    def static_alive(self, ids, time):
        """Mask of the static objects ids alive at time, False for an
        unknown id."""
        lives = [self._static.get(int(i), (0, 0, np.inf))[1:] for i in ids]
        return np.array([birth <= time < death for birth, death in lives], dtype=bool)

    def static_geometry(self, rows, time):
        """The static objects alive at time merged in one geometry, and the
        id of each of its cells.

        rows : (N, 8) id, x, y, z, q0, q1, q2, q3, as the static table
        without its time column.
        """
        rows = np.asarray(rows, dtype=np.float64).reshape(-1, 8)
        rows = rows[self.static_alive(rows[:, 0], time)]
        points, cells, ids = [], [], []
        n = 0
        for row in rows:
            geometry = self._static[int(row[0])][0]
            points.append(_rotate(row[4:8], geometry.points) + row[1:4])
            cells += [c + n for c in geometry.cells]
            ids.append(np.full(len(geometry.cells), int(row[0]), dtype=np.int64))
            n += len(geometry.points)
        points = np.vstack(points) if points else np.zeros((0, 3))
        ids = np.concatenate(ids) if ids else np.zeros(0, np.int64)
        return _Geometry(points, cells), ids

    def write(self, filename):
        """The glyph table, a MultiBlockDataSet with one PolyData per block."""
        version = np.array(VTKHDF_VERSION, dtype=np.int64)
        # the reader iterates the groups in creation order
        with h5py.File(filename, "w", track_order=True) as f:
            root = f.create_group("VTKHDF", track_order=True)
            root.attrs["Version"] = version
            root.attrs["Type"] = np.bytes_("MultiBlockDataSet")
            assembly = root.create_group("Assembly", track_order=True)
            for k, (name, geometry) in enumerate(self.blocks):
                block = f"Block_{k:06d}"
                g = root.create_group(block, track_order=True)
                g.attrs["Version"] = version
                g.attrs["Type"] = np.bytes_("PolyData")
                _write_polydata(g, geometry)
                # the assembly names must be unique
                assembly[f"{k}_{name}"] = h5py.SoftLink(f"/VTKHDF/{block}")


def _write_polydata(group, geometry):
    """A static PolyData: the cells of 2 points are lines, polygons else."""
    group.create_dataset("Points", data=geometry.points)
    group.create_dataset("NumberOfPoints", data=[len(geometry.points)])
    lines = [c for c in geometry.cells if len(c) == 2]
    polygons = [c for c in geometry.cells if len(c) > 2]
    for name in _TOPOLOGIES:
        cells = {"Lines": lines, "Polygons": polygons}.get(name, [])
        topology = group.create_group(name)
        connectivity = np.concatenate(cells) if cells else np.zeros(0, np.int64)
        offsets = np.concatenate(([0], np.cumsum([len(c) for c in cells])))
        topology.create_dataset("NumberOfCells", data=[len(cells)])
        topology.create_dataset("NumberOfConnectivityIds", data=[len(connectivity)])
        topology.create_dataset("Offsets", data=offsets.astype(np.int64))
        topology.create_dataset("Connectivity", data=connectivity.astype(np.int64))
    for name in ("PointData", "CellData", "FieldData"):
        group.create_group(name)


_PARAVIEW_STATE = """# ParaView state: paraview --state={state}
# or pvpython {state}, or File > Load State in ParaView.
# The bodies drawn with their shape: the 3D Glyphs representation picks in
# the glyph table (the shapes file) the block shape_index of each body,
# rotates it by the orientation quaternion and scales it by glyph_scale.
#
# Memory: the 3D Glyphs of ParaView up to 6.0 at least keep a copy of the
# glyph table for every step shown, never freed (fixed in VTK by commit
# e9652f91, "vtkOpenGLGlyph3DMapper: prevent excessive memory usage").
# Playing a long animation with many distinct shapes takes gigabytes:
# BODIES_AS_DISKS = True draws every body as a disk (a sphere in 3D) of
# radius disk_radius instead (in 2D, of the area of the shape), with the
# Point Gaussian representation, which keeps nothing.
BODIES_AS_DISKS = {disks}

from paraview.simple import *

bodies = OpenDataFile({bodies!r})
RenameSource("bodies", bodies)
view = GetActiveViewOrCreate("RenderView")
display = Show(bodies, view)
if BODIES_AS_DISKS:
    display.SetRepresentationType("Point Gaussian")
    display.ShaderPreset = {preset!r}
    display.ScaleByArray = 1
    display.SetScaleArray = ["POINTS", "disk_radius"]
    display.UseScaleFunction = 0
    display.GaussianRadius = 1.0
else:
    shapes = OpenDataFile({shapes!r})
    RenameSource("shapes", shapes)
    display.SetRepresentationType("3D Glyphs")
    display.GlyphType = shapes
    display.UseGlyphTable = 1
    display.UseCompositeGlyphTable = 1
    display.GlyphTableIndexArray = "shape_index"
    display.Orient = 1
    display.OrientationMode = "Quaternion"
    display.SelectOrientationVectors = "orientation"
    display.Scaling = 1
    display.ScaleMode = "Magnitude"
    display.SelectScaleArray = "glyph_scale"
    display.ScaleFactor = 1.0
ColorBy(display, ("POINTS", "velocity", "Magnitude"))
# the bodies start at rest: the range grows with the steps shown
velocity = GetColorTransferFunction("velocity")
velocity.AutomaticRescaleRangeMode = "Grow and update every timestep"
display.SetScalarBarVisibility(view, True)
GetAnimationScene().UpdateAnimationUsingDataTimeSteps()
{static}{view_2d}# the first Render resets the camera, it is set after
Render(view)
{camera}Render(view)
"""

_PARAVIEW_STATIC = """
# the static objects, their geometry merged, the id cell array tells
# which object a cell is of
static = OpenDataFile({static!r})
RenameSource("static", static)
static_display = Show(static, view)
static_display.SetRepresentationType("Surface")
static_display.ColorArrayName = ["POINTS", ""]
static_display.AmbientColor = static_display.DiffuseColor = [0.5, 0.5, 0.5]
static_display.LineWidth = 2.0
"""


_PARAVIEW_CONTACTS = """
# the contacts, hidden: one point per contact, colored by the reaction
# impulse; the normal and tangent parts are reaction_normal and
# reaction_tangent
contacts = OpenDataFile({contacts!r})
RenameSource("contacts", contacts)
contacts_display = Show(contacts, view)
contacts_display.SetRepresentationType("Point Gaussian")
contacts_display.GaussianRadius = {radius:.6g}
ColorBy(contacts_display, ("POINTS", "reaction", "Magnitude"))
reaction = GetColorTransferFunction("reaction")
reaction.AutomaticRescaleRangeMode = "Grow and update every timestep"
Hide(contacts, view)
"""


_PARAVIEW_FEM = """
# the finite element meshes, deformed, colored by the von Mises stress
fem = OpenDataFile({fem!r})
RenameSource("fem", fem)
fem_display = Show(fem, view)
fem_display.SetRepresentationType("Surface With Edges")
ColorBy(fem_display, ("CELLS", "von_mises"))
von_mises = GetColorTransferFunction("von_mises")
von_mises.AutomaticRescaleRangeMode = "Grow and update every timestep"
fem_display.SetScalarBarVisibility(view, True)
"""


def write_paraview_state(
    filename,
    bodies_filename,
    shapes_filename,
    dimension,
    static_filename=None,
    bounds=None,
    contacts_filename=None,
    typical_size=1.0,
    fem_filename=None,
    disks=False,
):
    """A ParaView python state drawing the bodies, and the static objects,
    with their shapes, and the contacts. disks: the default of
    BODIES_AS_DISKS in the state, the bodies drawn as disks."""
    # the glyph bounds are padded by the largest glyph: the camera is set on
    # the static objects and the bodies at their initial positions
    camera = "view.ResetCamera()\n"
    if bounds is not None:
        camera = f"view.ResetCamera({', '.join(f'{v:.6g}' for v in bounds)})\n"
    static = ""
    if static_filename is not None:
        static = _PARAVIEW_STATIC.format(static=os.path.abspath(static_filename))
    if fem_filename is not None:
        static += _PARAVIEW_FEM.format(fem=os.path.abspath(fem_filename))
    if contacts_filename is not None:
        static += _PARAVIEW_CONTACTS.format(
            contacts=os.path.abspath(contacts_filename),
            radius=0.25 * typical_size,
        )
    with open(filename, "w") as f:
        f.write(
            _PARAVIEW_STATE.format(
                state=os.path.basename(filename),
                bodies=os.path.abspath(bodies_filename),
                shapes=os.path.abspath(shapes_filename),
                static=static,
                disks=bool(disks),
                preset="Plain circle" if dimension == 2 else "Sphere",
                camera=camera,
                view_2d='view.InteractionMode = "2D"\n' if dimension == 2 else "",
            )
        )


Companions = namedtuple("Companions", "shapes state static contacts fem")


def companion_filenames(vtkhdf_filename):
    """The files written with a bodies VTKHDF file: the shapes, the ParaView
    state, the static objects, the contacts and the FEM meshes."""
    stem = os.path.splitext(vtkhdf_filename)[0]
    return Companions(
        f"{stem}_shapes.vtkhdf",
        f"{stem}_paraview.py",
        f"{stem}_static.vtkhdf",
        f"{stem}_contacts.vtkhdf",
        f"{stem}_fem.vtkhdf",
    )


def write_shapes_and_state(
    data, dimension, vtkhdf_filename, contacts=False, disks=False
):
    """The shapes file and the ParaView state of the bodies of a siconos
    file, returns the ShapeTable. contacts: draw the contacts file too.
    disks: the state draws the bodies as disks, for the large cases."""
    shapes = ShapeTable(data, dimension)
    if shapes.unsupported:
        print(
            "warning: shapes drawn as a cross, not supported: "
            + ", ".join(sorted(shapes.unsupported))
        )
    files = companion_filenames(vtkhdf_filename)
    shapes.fem = FemMeshes(data)
    if shapes.fem:
        shapes._extend_bounds(shapes.fem.points, 0.0)
    shapes.write(files.shapes)
    write_paraview_state(
        files.state,
        vtkhdf_filename,
        files.shapes,
        dimension,
        files.static if shapes.has_static else None,
        shapes.bounds,
        files.contacts if contacts else None,
        shapes.typical_size,
        files.fem if shapes.fem else None,
        disks,
    )
    return shapes


# --- finite element meshes ----------------------------------------------------
#
# The FEM bodies (input objects with a material) are drawn from their gmsh
# 2.2 mesh, kept in the ref group. Their displacements are stored in the
# siconos file in "spatial order", the mesh vertices sorted by (y, x), as
# mechanics_run builds it: the points of the meshes are put in this order.
# The strains and stresses are per T3 element, in the order of the mesh.


def _gmsh_triangles(text):
    """Nodes (N, 3) in file order and T3 triangles (M, 3) as indices in the
    nodes, from a gmsh 2.2 ASCII mesh."""
    lines = iter(text.splitlines())
    nodes, ids, triangles = [], [], []
    for line in lines:
        if line.startswith("$Nodes"):
            for _ in range(int(next(lines))):
                fields = next(lines).split()
                ids.append(int(fields[0]))
                nodes.append([float(v) for v in fields[1:4]])
        elif line.startswith("$Elements"):
            for _ in range(int(next(lines))):
                fields = [int(v) for v in next(lines).split()]
                # id, type, number of tags, tags..., nodes; type 2 = T3
                if fields[1] == 2:
                    triangles.append(fields[3 + fields[2] :])
    position = {node: k for k, node in enumerate(ids)}
    triangles = np.array([[position[n] for n in t] for t in triangles], np.int64)
    return np.array(nodes, dtype=np.float64), triangles.reshape(-1, 3)


class FemMeshes:
    """The meshes of the FEM bodies of a siconos mechanics HDF5 file, merged
    in one PolyData: points in spatial order, one body after the other."""

    def __init__(self, data):
        self.ids, self.points, self.triangles, self.body = [], [], [], []
        self.point_ranges, self.cell_ranges = {}, {}
        n_points = n_cells = 0
        for obj in data["input"].values():
            if "material" not in obj.attrs:
                continue
            contactors = [obj[c] for c in obj if "shape_name" in obj[c].attrs]
            mesh = data["ref"][contactors[0].attrs["shape_name"]][0]
            nodes, triangles = _gmsh_triangles(
                mesh.decode("utf-8") if isinstance(mesh, bytes) else str(mesh)
            )
            # the spatial order of mechanics_run, a stable sort on (y, x)
            order = np.lexsort((nodes[:, 0], nodes[:, 1]))
            rank = np.empty_like(order)
            rank[order] = np.arange(len(order))
            ds_id = int(obj.attrs["id"])
            self.ids.append(ds_id)
            self.points.append(nodes[order])
            self.triangles.append(rank[triangles] + n_points)
            self.body.append(np.full(len(triangles), ds_id, dtype=np.int64))
            self.point_ranges[ds_id] = slice(n_points, n_points + len(nodes))
            self.cell_ranges[ds_id] = slice(n_cells, n_cells + len(triangles))
            n_points += len(nodes)
            n_cells += len(triangles)
        if self.ids:
            self.points = np.vstack(self.points)
            self.triangles = np.vstack(self.triangles)
            self.body = np.concatenate(self.body)

    def __bool__(self):
        return len(self.ids) > 0

    def step_data(self, displacements, strains, stresses):
        """Points, point data and cell data of a step.

        displacements, strains, stresses : dicts ds id -> the row of the
        siconos table without its time and id columns (spatial order
        displacements, x and y interleaved; 3 components per element); a
        missing body has no displacement and NaN tensors.
        """
        n_points, n_cells = len(self.points), len(self.triangles)
        displacement = np.zeros((n_points, 3))
        strain = np.full((n_cells, 3), np.nan)
        stress = np.full((n_cells, 3), np.nan)
        for ds_id in self.ids:
            points, cells = self.point_ranges[ds_id], self.cell_ranges[ds_id]
            if ds_id in displacements:
                displacement[points, :2] = np.reshape(displacements[ds_id], (-1, 2))
            for table, values in ((strains, strain), (stresses, stress)):
                if ds_id in table:
                    values[cells] = np.reshape(table[ds_id], (-1, 3))
        sxx, syy, sxy = stress.T
        cell_data = {
            "body": self.body,
            "strain": strain,
            "stress": stress,
            # plane stress
            "von_mises": np.sqrt(sxx * sxx - sxx * syy + syy * syy + 3 * sxy * sxy),
        }
        return self.points + displacement, {"displacement": displacement}, cell_data


class FemWriter:
    """Append the deformed FEM meshes, step by step, to a VTKHDF file: a
    transient PolyData whose triangles are written once."""

    def __init__(self, filename, meshes, dtype=np.float32, append=False):
        self._dtype = dtype
        self._n_cells = len(meshes.triangles)
        if append and os.path.exists(filename):
            self._file = h5py.File(filename, "a")
            self._root = self._file["VTKHDF"]
            self._steps = self._root["Steps"]
            self._number_of_steps = int(self._steps.attrs["NSteps"])
            return
        self._file = h5py.File(filename, "w")
        self._number_of_steps = 0
        root = self._root = self._file.create_group("VTKHDF")
        root.attrs["Version"] = np.array(VTKHDF_VERSION, dtype=np.int64)
        root.attrs["Type"] = np.bytes_("PolyData")
        _extendable(root, "Points", (3,), dtype, 8192, None)
        _extendable(root, "NumberOfPoints", (), np.int64, _STEP_CHUNK_ROWS, None)
        for name in _TOPOLOGIES:
            topology = root.create_group(name)
            for count in ("NumberOfCells", "NumberOfConnectivityIds"):
                _extendable(topology, count, (), np.int64, _STEP_CHUNK_ROWS, None)
            if name == "Polygons":
                n = self._n_cells
                topology.create_dataset("Offsets", data=3 * np.arange(n + 1))
                topology.create_dataset(
                    "Connectivity", data=meshes.triangles.ravel().astype(np.int64)
                )
            else:
                topology.create_dataset("Offsets", data=np.zeros(1, np.int64))
                topology.create_dataset("Connectivity", data=np.zeros(0, np.int64))
        for name in ("PointData", "CellData", "FieldData"):
            root.create_group(name)
        steps = self._steps = root.create_group("Steps")
        steps.attrs["NSteps"] = 0
        _extendable(steps, "Values", (), np.float64, _STEP_CHUNK_ROWS, None)
        for name in ("PartOffsets", "NumberOfParts", "PointOffsets"):
            _extendable(steps, name, (), np.int64, _STEP_CHUNK_ROWS, None)
        for name in ("CellOffsets", "ConnectivityIdOffsets"):
            _extendable(
                steps, name, (len(_TOPOLOGIES),), np.int64, _STEP_CHUNK_ROWS, None
            )
        for name in ("PointDataOffsets", "CellDataOffsets", "FieldDataOffsets"):
            steps.create_group(name)

    def _array(self, group, name, values):
        """Append values to an array of the group, its step offset."""
        root, steps = self._root[group], self._steps[f"{group}Offsets"]
        if name not in root:
            dtype = np.int64 if np.issubdtype(values.dtype, np.integer) else self._dtype
            _extendable(root, name, values.shape[1:], dtype, 8192, None)
            _extendable(steps, name, (), np.int64, _STEP_CHUNK_ROWS, None)
        _append(steps[name], [_append(root[name], values)])

    def write_step(self, time, points, point_data, cell_data):
        root, steps = self._root, self._steps
        k = self._number_of_steps
        _append(steps["Values"], [time])
        _append(steps["PartOffsets"], [k])
        _append(steps["NumberOfParts"], [1])
        _append(steps["PointOffsets"], [_append(root["Points"], points)])
        _append(root["NumberOfPoints"], [len(points)])
        # every step reads the same triangles
        zeros = np.zeros((1, len(_TOPOLOGIES)), dtype=np.int64)
        _append(steps["CellOffsets"], zeros)
        _append(steps["ConnectivityIdOffsets"], zeros)
        for name in _TOPOLOGIES:
            cells = self._n_cells if name == "Polygons" else 0
            _append(root[name]["NumberOfCells"], [cells])
            _append(root[name]["NumberOfConnectivityIds"], [3 * cells])
        for name, values in point_data.items():
            self._array("PointData", name, np.asarray(values))
        for name, values in cell_data.items():
            self._array("CellData", name, np.asarray(values))
        self._number_of_steps += 1
        steps.attrs["NSteps"] = self._number_of_steps

    def flush(self):
        self._file.flush()

    def close(self):
        if self._file:
            self._file.close()
            self._file = None

    def __enter__(self):
        return self

    def __exit__(self, type_, value, traceback):
        self.close()


class StaticWriter:
    """Append the static objects, merged in one geometry, to a VTKHDF file:
    a transient PolyData of lines and polygons with the id of the object of
    each cell. A step is written only when the geometry changes: the
    reader shows the last step at or before the time asked.

    Unlike the bodies, the static objects are not glyphs: a glyph table of
    one block per static object costs ParaView memory at every step shown
    (see _PARAVIEW_STATE).
    """

    def __init__(self, filename, dtype=np.float32, compression=False, append=False):
        self._dtype = dtype
        self._compression = "gzip" if compression else None
        # the geometry of the last step written, unknown after a restart
        self._last = None
        if append and os.path.exists(filename):
            self._file = h5py.File(filename, "a")
            self._root = self._file["VTKHDF"]
            self._steps = self._root["Steps"]
            return
        self._file = h5py.File(filename, "w")
        root = self._root = self._file.create_group("VTKHDF")
        root.attrs["Version"] = np.array(VTKHDF_VERSION, dtype=np.int64)
        root.attrs["Type"] = np.bytes_("PolyData")
        self._extendable(root, "Points", (3,), dtype)
        self._per_step(root, "NumberOfPoints", (), np.int64)
        for name in _TOPOLOGIES:
            topology = root.create_group(name)
            for count in ("NumberOfCells", "NumberOfConnectivityIds"):
                self._per_step(topology, count, (), np.int64)
            # a step without cell reads the first offset
            _append(self._extendable(topology, "Offsets", (), np.int64), [0])
            self._extendable(topology, "Connectivity", (), np.int64)
        root.create_group("PointData")
        self._extendable(root.create_group("CellData"), "id", (), np.int64)
        root.create_group("FieldData")
        steps = self._steps = root.create_group("Steps")
        steps.attrs["NSteps"] = 0
        self._per_step(steps, "Values", (), np.float64)
        for name in ("PartOffsets", "NumberOfParts", "PointOffsets"):
            self._per_step(steps, name, (), np.int64)
        for name in ("CellOffsets", "ConnectivityIdOffsets"):
            self._per_step(steps, name, (len(_TOPOLOGIES),), np.int64)
        steps.create_group("PointDataOffsets")
        self._per_step(steps.create_group("CellDataOffsets"), "id", (), np.int64)
        steps.create_group("FieldDataOffsets")

    def _extendable(self, group, name, shape, dtype):
        return _extendable(group, name, shape, dtype, 8192, self._compression)

    def _per_step(self, group, name, shape, dtype):
        return _extendable(group, name, shape, dtype, _STEP_CHUNK_ROWS, None)

    def _cell_offset(self, topology):
        """The CellOffsets of a new step and the row of Offsets where its
        offsets go, padding Offsets up to it.

        vtkHDFReader reads the offsets of a step at CellOffsets + the
        number of the previous steps whose cells start before
        CellOffsets, not at the true row when a step without cell comes
        just before: the step is written at the first row it can read.
        """
        offsets = topology["Offsets"]
        counts = topology["NumberOfCells"][:]
        # the first cell of each previous step
        starts = np.concatenate(([0], np.cumsum(counts)[:-1])) if len(counts) else []
        row = offsets.shape[0]

        def read_row(cell_offset):
            return cell_offset + int(np.searchsorted(starts, cell_offset, "left"))

        low, high = 0, row
        while low < high:
            middle = (low + high) // 2
            if read_row(middle) >= row:
                high = middle
            else:
                low = middle + 1
        padding = read_row(low) - row
        if padding:
            _append(offsets, np.zeros(padding, np.int64))
        return low

    def write_step(self, time, geometry, ids):
        """Append one time step, unless the geometry is the one of the
        last step written.

        geometry : _Geometry, its cells of 2 points are lines, polygons
            else.
        ids : (number of cells,) the id of the object of each cell.
        """
        lines = [c for c in geometry.cells if len(c) == 2]
        polygons = [c for c in geometry.cells if len(c) > 2]
        is_line = np.array([len(c) == 2 for c in geometry.cells], dtype=bool)
        ids = np.asarray(ids, dtype=np.int64)
        # the cell data in the order of the topologies
        ids = np.concatenate((ids[is_line], ids[~is_line]))
        cells = {"Lines": lines, "Polygons": polygons}
        connectivity = {
            name: np.concatenate(c) if c else np.zeros(0, np.int64)
            for name, c in cells.items()
        }
        current = (
            geometry.points,
            connectivity["Lines"],
            connectivity["Polygons"],
            ids,
        )
        if self._last is not None and all(
            np.array_equal(a, b) for a, b in zip(current, self._last)
        ):
            return
        self._last = current

        root, steps = self._root, self._steps
        k = int(steps.attrs["NSteps"])
        _append(steps["Values"], [time])
        _append(steps["PartOffsets"], [k])
        _append(steps["NumberOfParts"], [1])
        points = geometry.points.astype(self._dtype, copy=False)
        _append(steps["PointOffsets"], [_append(root["Points"], points)])
        _append(root["NumberOfPoints"], [len(points)])
        cell_offsets, connectivity_offsets = [], []
        for name in _TOPOLOGIES:
            topology = root[name]
            step_cells = cells.get(name, [])
            step_connectivity = connectivity.get(name, np.zeros(0, np.int64))
            if step_cells:
                cell_offsets.append(self._cell_offset(topology))
                counts = np.cumsum([len(c) for c in step_cells])
                _append(topology["Offsets"], np.concatenate(([0], counts)))
            else:
                cell_offsets.append(0)
            connectivity_offsets.append(
                _append(topology["Connectivity"], step_connectivity)
            )
            _append(topology["NumberOfCells"], [len(step_cells)])
            _append(topology["NumberOfConnectivityIds"], [len(step_connectivity)])
        _append(steps["CellOffsets"], [cell_offsets])
        _append(steps["ConnectivityIdOffsets"], [connectivity_offsets])
        offset = _append(root["CellData"]["id"], ids)
        _append(steps["CellDataOffsets"]["id"], [offset])
        steps.attrs["NSteps"] = k + 1

    def flush(self):
        self._file.flush()

    def close(self):
        if self._file:
            self._file.close()
            self._file = None

    def __enter__(self):
        return self

    def __exit__(self, type_, value, traceback):
        self.close()


def fem_rows(data, time=None):
    """The rows of the FEM tables, dicts ds id -> row without the time and
    id columns: displacements, strains, stresses. At time, or the last
    rows if time is None; a table without a row at time is left out."""
    tables = ({}, {}, {})
    for name in data:
        for prefix, table in zip(
            ("fem_displacements_", "fem_epsilon_", "fem_sigma_"), tables
        ):
            if name.startswith(prefix) and data[name].shape[0] > 0:
                row = data[name][-1]
                if time is None or abs(row[0] - time) <= 1e-9 * max(1, abs(time)):
                    table[int(name[len(prefix) :])] = row[2:]
    return tables


def bodies_point_data(positions, velocities, radii=None):
    """Points and point arrays of the bodies from the siconos IO arrays.

    positions : (N, 1 + 3) in 2D (id, x, y, theta) or (N, 1 + 7) in 3D
        (id, x, y, z, q0, q1, q2, q3), as io.positions(nsds).
    velocities : (N, 1 + 3) in 2D (id, vx, vy, w) or (N, 1 + 6) in 3D,
        as io.velocities(nsds), possibly with extra trailing columns.
    radii : (N, 2) (id, r), as io.radii(nsds), optional.

    The rows of the three arrays must be in the same body order.
    """
    positions = np.asarray(positions)
    velocities = np.asarray(velocities)
    n = positions.shape[0]
    ids = positions[:, 0]

    if positions.shape[1] == 4:
        points = np.column_stack((positions[:, 1:3], np.zeros(n)))
        half = positions[:, 3] / 2.0
        orientation = np.zeros((n, 4))
        orientation[:, 0] = np.cos(half)
        orientation[:, 3] = np.sin(half)
        velocity = np.column_stack((velocities[:, 1:3], np.zeros(n)))
        angular_velocity = np.zeros((n, 3))
        angular_velocity[:, 2] = velocities[:, 3]
    else:
        points = positions[:, 1:4]
        orientation = positions[:, 4:8]
        velocity = velocities[:, 1:4]
        angular_velocity = velocities[:, 4:7]

    point_data = {
        "id": ids,
        "orientation": orientation,
        "velocity": velocity,
        "angular_velocity": angular_velocity,
    }
    if radii is not None:
        radii = np.asarray(radii)
        if radii.shape[0] != n or not np.array_equal(radii[:, 0], ids):
            raise ValueError("radii and positions are not in the same body order")
        point_data["radius"] = radii[:, 1]
    return points, point_data


# the cf table columns, after its time column
_CF = {
    "mu": 0,
    "point_a": slice(1, 4),
    "point_b": slice(4, 7),
    "normal": slice(7, 10),
    "reaction": slice(10, 13),
    "gap": 13,
    "relative_velocity": slice(16, 19),
    "reaction_local": slice(19, 22),
    "interaction": 22,
    "ds1": 23,
    "ds2": 24,
}


def contacts_point_data(rows):
    """Points and point arrays of the contacts.

    rows : (N, 25) the rows of the cf table without the time column. A
    contact is drawn at the middle of its points A and B.
    """
    rows = np.asarray(rows, dtype=np.float64).reshape(-1, 25)
    points = 0.5 * (rows[:, _CF["point_a"]] + rows[:, _CF["point_b"]])
    local = rows[:, _CF["reaction_local"]]
    point_data = {
        "normal": rows[:, _CF["normal"]],
        "reaction": rows[:, _CF["reaction"]],
        "reaction_normal": local[:, 0],
        "reaction_tangent": np.linalg.norm(local[:, 1:], axis=1),
        "gap": rows[:, _CF["gap"]],
        "relative_velocity": rows[:, _CF["relative_velocity"]],
        "mu": rows[:, _CF["mu"]],
    }
    for name in ("interaction", "ds1", "ds2"):
        point_data[name] = rows[:, _CF[name]].astype(np.int64)
    return points, point_data


def siconos_rows(time, points, point_data):
    """The rows of the siconos dynamic and velocities tables of a step.

    The inverse of bodies_point_data on the tables: (N, 9) rows time,
    id, x, y, z, q0, q1, q2, q3 and (N, 8) rows time, id, v, w.
    """
    n = points.shape[0]
    times = np.full((n, 1), time)
    ids = np.asarray(point_data["id"], dtype=np.float64).reshape(n, 1)
    dynamic = np.hstack((times, ids, points, point_data["orientation"]))
    velocities = np.hstack(
        (times, ids, point_data["velocity"], point_data["angular_velocity"])
    )
    return dynamic, velocities


def _step_bounds(table, block_rows):
    """(time, first row, end row) of the steps of a siconos table.

    The rows of a step are contiguous. The time column is read by
    blocks, never as a whole.
    """
    n = table.shape[0]
    starts, times = [], []
    previous = None
    for first in range(0, n, block_rows):
        block = table[first : first + block_rows, 0]
        changes = np.flatnonzero(block[1:] != block[:-1]) + 1
        if previous is None or block[0] != previous:
            changes = np.concatenate(([0], changes))
        starts.extend(first + changes)
        times.extend(block[changes])
        previous = block[-1]
    ends = [*starts[1:], n]
    return [(float(t), int(s), int(e)) for t, s, e in zip(times, starts, ends)]


def _indexed_steps(data):
    """The steps of a file with a step index.

    Yields (time, dynamic rows, velocities rows, radii rows or None) for
    the steps with bodies. A step without radii (the initial output)
    takes those of the next step that has some.
    """
    index = data["step_index"]
    dynamic = index["dynamic"][:].astype(np.int64)
    times = index["dynamic"][:, 0]
    velocities = index["velocities"][:].astype(np.int64)
    radii = index["radii"][:].astype(np.int64) if "radii" in data else None
    with_radii = radii is not None and radii[:, 2].sum() > 0

    def rows(table, k):
        return slice(table[k, 1], table[k, 1] + table[k, 2])

    for k in np.flatnonzero(dynamic[:, 2] > 0):
        step_radii = None
        if with_radii:
            later = np.flatnonzero(radii[k:, 2] > 0)
            nearest = k + later[0] if later.size else np.flatnonzero(radii[:, 2])[-1]
            step_radii = rows(radii, nearest)
        yield float(times[k]), rows(dynamic, k), rows(velocities, k), step_radii


def _time_column_steps(data, block_rows):
    """The steps of a file without step index, from the time column.

    The radii table, which has no time column and no entry for the
    initial output, is used only when its size lines up with dynamic.
    """
    dynamic = data["dynamic"]
    bounds = _step_bounds(dynamic, block_rows)
    radii = data.get("radii")
    shift = bounds[0][2] if bounds else 0
    if radii is not None and (
        radii.shape[0] == 0 or radii.shape[0] != dynamic.shape[0] - shift
    ):
        print(
            f"warning: {radii.shape[0]} radii rows do not match "
            f"{dynamic.shape[0]} dynamic rows, radii not written"
        )
        radii = None
    for k, (time, start, end) in enumerate(bounds):
        step_radii = None
        if radii is not None:
            # the initial output takes the radii of the next one
            _, first, last = bounds[min(max(k, 1), len(bounds) - 1)]
            step_radii = slice(first - shift, last - shift)
        yield time, slice(start, end), slice(start, end), step_radii


def convert(
    io_filename, vtkhdf_filename, compression=False, block_rows=1 << 22, disks=False
):
    """Write the bodies of a siconos mechanics HDF5 file in a VTKHDF file.

    The dynamic and velocities tables are read one step at a time, with
    the step index of the file when it has one. The radii table (vnative
    backend) is used when it lines up with them. The shapes file and the
    ParaView state are written alongside, the state draws the bodies as
    disks if disks.
    Returns the number of steps written.
    """
    with h5py.File(io_filename, "r") as io:
        data = io["data"]
        dynamic = data["dynamic"]
        velocities = data["velocities"]
        if velocities.shape[0] != dynamic.shape[0]:
            raise ValueError(
                f"{velocities.shape[0]} velocities rows for {dynamic.shape[0]} dynamic rows"
            )
        if "step_index" in data:
            steps = _indexed_steps(data)
        else:
            steps = _time_column_steps(data, block_rows)
        dimension = int(io.attrs.get("dimension", 3))
        cf = data.get("cf")
        with_contacts = cf is not None and cf.shape[0] > 0
        shapes = write_shapes_and_state(
            data, dimension, vtkhdf_filename, contacts=with_contacts, disks=disks
        )
        files = companion_filenames(vtkhdf_filename)

        number_of_steps = 0
        with BodiesWriter(vtkhdf_filename, compression=compression) as writer:
            for time, dynamic_rows, velocities_rows, radii_rows in steps:
                positions = dynamic[dynamic_rows, 1:]
                step_radii = None
                if radii_rows is not None:
                    step_radii = data["radii"][radii_rows]
                    if not np.array_equal(step_radii[:, 0], positions[:, 0]):
                        step_radii = _radii_of(positions[:, 0], step_radii)
                points, point_data = bodies_point_data(
                    positions, velocities[velocities_rows, 1:], step_radii
                )
                point_data.update(shapes.point_data(point_data["id"]))
                writer.write_step(time, points, point_data)
                number_of_steps += 1

        if with_contacts:
            _write_contacts(data, files.contacts, compression, block_rows)
        if shapes.fem:
            _write_fem(data, files.fem, shapes.fem)

        if shapes.has_static:
            static_filename = files.static
            static = data["static"]
            if "step_index" in data:
                rows = data["step_index"]["static"][:]
                bounds = [(float(t), int(s), int(s + n)) for t, s, n in rows if n > 0]
            else:
                bounds = _step_bounds(static, block_rows)
            with StaticWriter(static_filename, compression=compression) as writer:
                for time, start, end in bounds:
                    writer.write_step(
                        time, *shapes.static_geometry(static[start:end, 1:], time)
                    )
    return number_of_steps


def _write_contacts(data, filename, compression, block_rows):
    """The contacts of the cf table, one step per bodies step, empty or not:
    ParaView would keep showing the contacts of the last step written."""
    cf = data["cf"]
    if "step_index" in data:
        index = data["step_index"]
        bodies = index["dynamic"][:]
        rows = index["cf"][:].astype(np.int64)
        steps = [
            (float(t), slice(rows[k, 1], rows[k, 1] + rows[k, 2]))
            for k, t in enumerate(bodies[:, 0])
            if bodies[k, 2] > 0
        ]
    else:
        cf_rows = {t: slice(s, e) for t, s, e in _step_bounds(cf, block_rows)}
        times = [t for t, _, _ in _step_bounds(data["dynamic"], block_rows)]
        steps = [(t, cf_rows.get(t, slice(0, 0))) for t in times]
    with BodiesWriter(filename, compression=compression) as writer:
        for time, rows in steps:
            writer.write_step(time, *contacts_point_data(cf[rows, 1:]))


def _write_fem(data, filename, meshes):
    """The FEM meshes at each row of the FEM tables, after an undeformed
    step at the first output time if the tables start later."""
    tables = [
        {i: data[f"{prefix}{i}"] for i in meshes.ids if f"{prefix}{i}" in data}
        for prefix in ("fem_displacements_", "fem_epsilon_", "fem_sigma_")
    ]
    rows = max((t.shape[0] for table in tables for t in table.values()), default=0)
    fem_first = min(
        (float(t[0, 0]) for t in tables[0].values() if len(t)), default=None
    )
    first = float(data["dynamic"][0, 0]) if data["dynamic"].shape[0] else None
    with FemWriter(filename, meshes) as writer:
        if first is not None and (fem_first is None or first < fem_first):
            writer.write_step(first, *initial_fem_step(meshes))
        for k in range(rows):
            step = [
                {i: t[k, 2:] for i, t in table.items() if k < t.shape[0]}
                for table in tables
            ]
            time = next(float(t[k, 0]) for t in tables[0].values() if k < t.shape[0])
            writer.write_step(time, *meshes.step_data(*step))


def initial_fem_step(meshes):
    """The undeformed meshes, no strain nor stress."""
    zero = {i: np.zeros(3 * (r.stop - r.start)) for i, r in meshes.cell_ranges.items()}
    return meshes.step_data({}, zero, zero)


def _radii_of(ids, radii):
    """(id, r) rows for the bodies ids, r = 0 for a body not in radii."""
    order = np.argsort(radii[:, 0])
    known = radii[order, 0]
    where = np.clip(np.searchsorted(known, ids), 0, max(len(known) - 1, 0))
    found = known[where] == ids
    return np.column_stack((ids, np.where(found, radii[order, 1][where], 0.0)))
