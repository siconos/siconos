# Siconos is a program dedicated to modeling, simulation and control
# of non smooth dynamical systems.
#
# Copyright 2025 INRIA.
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

"""Functions and tools used in mechanics_run and mechanics_hdf5"""

import tempfile
import sys
import shutil
import subprocess
import os
from contextlib import contextmanager
import numpy as np
import time
import warnings
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
    (_fid, tfilename) = tempfile.mkstemp(suffix=suffix, prefix=prefix)
    fid = open(tfilename, "w")
    if contents is not None:
        fid.write(contents)
        fid.flush()

    class TmpFile:

        def __init__(self, fid, name):
            self.fid = fid
            self.name = name

        def __getitem__(self, n):
            if n == 0:
                return self.fid
            elif n == 1:
                return self.name
            else:
                raise IndexError

    r = TmpFile(fid, tfilename)

    yield r
    fid.close()
    os.close(_fid)
    if not debug:
        os.remove(tfilename)


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
        for i in range(0, points.GetNumberOfTuples()):
            p = points.GetTuple(i)
            apoints[0, i] = p[0]
            apoints[1, i] = p[1]
            apoints[2, i] = p[2]

        if scale is not None:
            apoints *= scale

        aindices = np.empty(num_triangles * 3, dtype=int)

        for i in range(0, num_triangles):
            c = polydata.GetCell(i)
            aindices[i * 3 + 0] = c.GetPointIds().GetId(0)
            aindices[i * 3 + 1] = c.GetPointIds().GetId(1)
            aindices[i * 3 + 2] = c.GetPointIds().GetId(2)

        shape = siconos.mechanics.collision.SiconosMesh(list(aindices), apoints)
        dims = apoints.max(axis=0) - apoints.min(axis=0)

    else:  # assume convex shape
        coors = dict()
        for i in range(0, points.GetNumberOfTuples()):
            coors[points.GetTuple(i)] = 1
        coors = np.array(coors.keys())
        dims = coors.max(axis=0) - coors.min(axis=0)
        shape = siconos.mechanics.collision.SiconosConvexHull(coors.keys())

    return shape, dims

def extract_bc_global_dofs(fesolid, mesh_data):
    """
    Return global DOF indices for vertices with physical tag 2 ("Dirichlet BC").
    """
    import meshio
    import tempfile

    with tempfile.NamedTemporaryFile(suffix='.msh', mode='w') as tmp:
        tmp.write(mesh_data)
        tmp.flush()
        mesh = meshio.read(tmp.name)

    # 1. "Dirichlet BC"
    tag_to_name = {int(tag): name.decode() if isinstance(name, bytes) else name 
                   for name, (tag, dim) in mesh.field_data.items()}
    dirichlet_tag = next((tag for tag, name in tag_to_name.items() if name == "Dirichlet BC"), None)
    if dirichlet_tag is None:
        return []

    # 2. collect boundary vertex IDs (1-based GMSH indices) with tag 2
    bc_vertex_ids = set()
    for cell_block, data_block in zip(mesh.cells, mesh.cell_data.get('gmsh:physical', [])):
        if cell_block.type == 'vertex':
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
        if cell_block.type != 'triangle':
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
