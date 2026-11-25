# Siconos is a program dedicated to modeling, simulation and control
# of non smooth dynamical systems.
#
# SPDX-FileCopyrightText: INRIA
# SPDX-License-Identifier: Apache-2.0
#

"""Functions and tools used in mechanics_run and mechanics_hdf5"""

import os
import shutil
import subprocess
import sys
import tempfile
import time
import warnings
from contextlib import contextmanager

import numpy as np
import siconos.mechanics.collision

warnings.simplefilter("always", UserWarning)


def arguments():
    """Returns tuple containing dictionary of calling function's
    named arguments and a list of calling function's unnamed
    positional arguments.
    """
    from inspect import getargvalues, stack

    posname, kwname, args = getargvalues(stack()[1][0])[-3:]
    posargs = args.pop(posname, [])
    args.update(args.pop(kwname, []))
    return args, posargs


@contextmanager
def tmpfile(suffix="", prefix="siconos_io", contents=None, debug=False):
    """
    A context manager for a named temporary file.
    """
    fd, filename = tempfile.mkstemp(
        suffix=suffix,
        prefix=prefix,
    )
    fid = os.fdopen(fd, "w")

    try:
        if contents is not None:
            fid.write(contents)
            fid.flush()

        class TmpFile:
            def __init__(self, fid, name):
                self.fid = fid
                self.name = name

        yield TmpFile(fid, filename)

    finally:
        fid.close()
        if not debug:
            os.remove(filename)


def warn(msg):
    warnings.warn(msg, category=UserWarning)


def object_id(obj):
    """returns an unique object identifier"""
    return obj.__hash__()


def str_of_file(filename):
    with open(filename, "r") as f:
        return str(f.read())


def file_of_str(filename, string):
    if not os.path.exists(os.path.dirname(filename)):
        try:
            os.makedirs(os.path.dirname(filename))
        except OSError as exc:
            if exc.errno != exc.errno.EEXIST:
                raise

    with open(filename, "w") as f:
        f.write(string)


def get_open_fds() -> int:
    """Get the number of open file descriptors for the current process."""
    lsof_path = shutil.which("lsof")
    if lsof_path is None:
        raise NotImplementedError("Didn't handle unavailable lsof.")
    raw_procs = subprocess.check_output(
        [lsof_path, "-w", "-Ff", "-p", str(os.getpid())]
    )

    def filter_fds(lsof_entry: str) -> bool:
        return lsof_entry.startswith("f") and lsof_entry[1:].isdigit()

    fds = list(filter(filter_fds, raw_procs.decode().split(os.linesep)))

    return len(fds)


time_measure = time.perf_counter
if sys.version_info.major + 0.1 * sys.version_info.minor < 3.3:
    time_measure = time.clock


class Timer:
    def __init__(self):
        self._t0 = time_measure()

    def elapsed(self):
        return time_measure() - self._t0

    def update(self):
        self._t0 = time_measure()


#
# load .vtp file
#
def load_siconos_mesh(shape_filename, scale=None):
    """
    loads a vtk .vtp file and returns a SiconosMesh shape
    WARNING triangles cells assumed!
    """
    import vtk

    reader = vtk.vtkXMLPolyDataReader()
    reader.SetFileName(shape_filename)
    reader.Update()

    polydata = reader.GetOutput()
    points = polydata.GetPoints().GetData()
    num_points = points.GetNumberOfTuples()
    num_triangles = polydata.GetNumberOfCells()

    shape = None

    if polydata.GetCellType(0) == 5:
        apoints = np.empty((3, num_points), dtype=np.float64, order="F")
        for i in range(points.GetNumberOfTuples()):
            p = points.GetTuple(i)
            apoints[0, i] = p[0]
            apoints[1, i] = p[1]
            apoints[2, i] = p[2]

        if scale is not None:
            apoints *= scale

        aindices = np.empty(num_triangles * 3, dtype=int)

        for i in range(num_triangles):
            c = polydata.GetCell(i)
            aindices[i * 3 + 0] = c.GetPointIds().GetId(0)
            aindices[i * 3 + 1] = c.GetPointIds().GetId(1)
            aindices[i * 3 + 2] = c.GetPointIds().GetId(2)

        shape = siconos.mechanics.collision.SiconosMesh(list(aindices), apoints)
        dims = apoints.max(axis=0) - apoints.min(axis=0)

    else:  # assume convex shape
        coors = {}
        for i in range(points.GetNumberOfTuples()):
            coors[points.GetTuple(i)] = 1
        coors = np.array(coors.keys())
        dims = coors.max(axis=0) - coors.min(axis=0)
        shape = siconos.mechanics.collision.SiconosConvexHull(coors.keys())

    return shape, dims


def extract_bc_global_dofs(fesolid, mesh_data):
    """
    Return global DOF indices for vertices with physical tag 2 ("Dirichlet BC").
    """
    import tempfile

    import meshio

    with tempfile.NamedTemporaryFile(suffix=".msh", mode="w") as tmp:
        tmp.write(mesh_data)
        tmp.flush()
        mesh = meshio.read(tmp.name)

    # 1. "Dirichlet BC"
    tag_to_name = {
        int(tag): name.decode() if isinstance(name, bytes) else name
        for name, (tag, dim) in mesh.field_data.items()
    }
    dirichlet_tag = next(
        (tag for tag, name in tag_to_name.items() if name == "Dirichlet BC"), None
    )
    if dirichlet_tag is None:
        return []

    # 2. collect boundary vertex IDs (1-based GMSH indices) with tag 2
    bc_vertex_ids = set()
    for cell_block, data_block in zip(
        mesh.cells, mesh.cell_data.get("gmsh:physical", [])
    ):
        if cell_block.type == "vertex":
            tags = data_block[0] if isinstance(data_block, list) else data_block
            for i, tag in enumerate(tags):
                if tag == dirichlet_tag:
                    bc_vertex_ids.add(cell_block.data[i][0])  # 1-based

    if not bc_vertex_ids:
        return []

    # 3. simulate FEModel::init() registration order:
    #    traverse elements in file order; for each element's vertices,
    #    if not yet registered, assign next DOF block.
    vertex_to_dof = {}  # 1-based vertex_id -> (dof_x, dof_y)
    dof_counter = 0
    dim = 2  # 2D FEM

    for cell_block in mesh.cells:
        if cell_block.type != "triangle":
            continue
        for tri in cell_block.data:
            for vid_1based in tri:
                if vid_1based not in vertex_to_dof:
                    vertex_to_dof[vid_1based] = (dof_counter, dof_counter + 1)
                    dof_counter += dim

    # 4. extract global DOFs for boundary vertices in sorted order
    #    (sorted by vertex_id for determinism)
    global_dofs = []
    for vid in sorted(bc_vertex_ids):
        if vid in vertex_to_dof:
            global_dofs.extend(vertex_to_dof[vid])

    return global_dofs


def chain_contact_segments(segments):
    """Order the contact segments of a FEM mesh into one chain.

    segments: the contact line elements, in any order and orientation, each
    as the dofs (x0, y0, x1, y1) of its two vertices
    (FiniteElementModel.getContactSegments).

    Returns (vertices, chained): vertices, the x dof of each vertex along the
    chain, the first one repeated at the end if the chain is closed;
    chained, the segments in chain order, oriented so that segment i goes
    from vertices[i] to vertices[i + 1]. This is the layout of the vnative
    chained_segment shape: n vertices for n - 1 segments, the dofs of the
    segment i at 4 i in the global indices.

    Raises ValueError if the segments do not form a single chain (a branch,
    or several disjoint chains)."""

    segments = [tuple(int(d) for d in s) for s in segments]
    if not segments:
        return [], []

    # the segments at each vertex, a vertex being known by its x dof
    at_vertex = {}
    for k, s in enumerate(segments):
        for v in (s[0], s[2]):
            at_vertex.setdefault(v, []).append(k)

    branches = [v for v, ks in at_vertex.items() if len(ks) > 2]
    if branches:
        raise ValueError(
            f"contact segments: {len(branches)} vertices shared by more than two "
            "segments, the contact boundary is not a single chain"
        )

    # an open chain starts at one of its two ends, a closed one anywhere
    ends = [v for v, ks in at_vertex.items() if len(ks) == 1]
    start = ends[0] if ends else segments[0][0]

    vertices = [start]
    chained = []
    used = set()
    v = start
    while True:
        free = [k for k in at_vertex[v] if k not in used]
        if not free:
            break
        k = free[0]
        used.add(k)
        s = segments[k]
        if s[0] != v:
            s = (s[2], s[3], s[0], s[1])
        chained.append(s)
        v = s[2]
        vertices.append(v)

    if len(chained) != len(segments):
        raise ValueError(
            f"contact segments: {len(segments) - len(chained)} of {len(segments)} "
            "segments are not connected to the chain, the contact boundary must be "
            "a single closed loop or open line"
        )

    return vertices, chained
