#!/usr/bin/env @Python_EXECUTABLE@
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

"""Tools for mechanics_run or hdf5, based on occ API."""

import os
import numpy as np

import vtk

from OCC.Core.GProp import GProp_GProps
from OCC.Core.gp import gp_Ax1, gp_Dir
import siconos.mechanics.occ
from OCC.Core.TopAbs import TopAbs_FACE
from OCC.Core.TopAbs import TopAbs_EDGE
from OCC.Core.TopExp import TopExp_Explorer
from OCC.Core.TopoDS import topods_Face, topods_Edge
from OCC.Core.STEPControl import STEPControl_Reader
from OCC.Core.IGESControl import IGESControl_Reader
from OCC.Core.BRep import BRep_Builder
from OCC.Core.TopoDS import TopoDS_Compound
from OCC.Core.IFSelect import IFSelect_RetDone, IFSelect_ItemsByEntity
from OCC.Core.StlAPI import StlAPI_Writer
from OCC.Core.BRepMesh import BRepMesh_IncrementalMesh
from OCC.Core.BRepTools import BRepTools_ShapeSet
import OCC.Core.BRepGProp

import siconos.io.tools


#
# inertia
#
def compute_inertia_and_center_of_mass(shapes, io=None):
    """
    Computes inertia from a list of volumes

    Returns
    -------
    mass
    center_of_mass
    inertia
    inertia_matrix
    """

    system = GProp_GProps()

    for shape in shapes:

        iprops = GProp_GProps()

        if shape.data is None:
            if io is not None:
                shape.data = io._shape.get(shape.shape_name, new_instance=True)
            else:
                siconos.io.tools.warn("cannot get shape {0}".format(shape.shape_name))
                return None

        ishape = siconos.mechanics.occ.OccContactShape(shape.data)
        siconos.mechanics.occ.occ_move(
            ishape, np.concatenate([shape.translation, shape.orientation], axis=0)
        )

        OCC.Core.BRepGProp.brepgprop.VolumeProperties(shape.data, iprops)

        density = None

        if hasattr(shape, "mass") and shape.mass is not None:
            density = shape.mass / iprops.Mass()

        elif shape.parameters is not None and hasattr(shape.parameters, "density"):
            density = shape.parameters.density
            # print('shape.parameters.density:', shape.parameters.density)
        else:
            density = 1.0

        assert density is not None
        # print("shape", shape.shape_name)
        # print('density:', density)
        # print('iprops.Mass():', iprops.Mass())

        system.Add(iprops, density)

    mass = system.Mass()
    assert system.Mass() > 0.0

    computed_com = system.CentreOfMass()

    gp_mat = system.MatrixOfInertia()
    inertia_matrix = np.zeros((3, 3), dtype=np.float64)
    for i in range(0, 3):
        for j in range(0, 3):
            inertia_matrix[i, j] = gp_mat.Value(i + 1, j + 1)

    I1 = system.MomentOfInertia(gp_Ax1(computed_com, gp_Dir(1, 0, 0)))
    I2 = system.MomentOfInertia(gp_Ax1(computed_com, gp_Dir(0, 1, 0)))
    I3 = system.MomentOfInertia(gp_Ax1(computed_com, gp_Dir(0, 0, 1)))

    inertia = [I1, I2, I3]
    center_of_mass = np.array(
        [computed_com.Coord(1), computed_com.Coord(2), computed_com.Coord(3)],
        dtype=np.float64,
    )

    return mass, center_of_mass, inertia, inertia_matrix


def occ_topo_list(shape):
    """return the edges & faces from `shape`

    :param shape: a TopoDS_Shape
    :return: a list of edges and faces
    """

    topExp = TopExp_Explorer()
    topExp.Init(shape, TopAbs_FACE)
    faces = []
    edges = []

    while topExp.More():
        face = topods_Face(topExp.Current())
        faces.append(face)
        topExp.Next()

    topExp.Init(shape, TopAbs_EDGE)

    while topExp.More():
        edge = topods_Edge(topExp.Current())
        edges.append(edge)
        topExp.Next()

    return faces, edges


def occ_load_file(filename):
    """
    Build a TopoDS from a file (step or igs)


    Parameters
    ----------
    filename : input file name (must be "stp", "step", "igs", or "iges"

    Returns
    -------
    TopoDS_Compound
        Compound shape built from all the shapes found in the IGES file.
    """
    reader_switch = {
        "stp": STEPControl_Reader,
        "step": STEPControl_Reader,
        "igs": IGESControl_Reader,
        "iges": IGESControl_Reader,
    }

    builder = BRep_Builder()
    comp = TopoDS_Compound()
    builder.MakeCompound(comp)
    reader = reader_switch[os.path.splitext(filename)[1][1:].lower()]()

    status = reader.ReadFile(filename)
    if status != IFSelect_RetDone:
        raise RuntimeError("occ_load_file read failed")

    failsonly = False
    reader.PrintCheckLoad(failsonly, IFSelect_ItemsByEntity)
    reader.PrintCheckTransfer(failsonly, IFSelect_ItemsByEntity)

    reader.TransferRoots()
    nbs = reader.NbShapes()

    for i in range(1, nbs + 1):
        shape = reader.Shape(i)
        builder.Add(comp, shape)

    return comp


def topods_shape_reader(shape, deflection=0.001):

    stl_writer = StlAPI_Writer()
    reader = vtk.vtkSTLReader()

    with siconos.io.tools.tmpfile(suffix=".stl") as tmpf:
        mesh = BRepMesh_IncrementalMesh(shape, deflection)
        mesh.Perform()
        assert mesh.IsDone()
        stl_writer.SetASCIIMode(False)
        stl_writer.Write(shape, tmpf[1])
        tmpf[0].flush()

        reader.SetFileName(tmpf[1])
        reader.Update()

    return reader


def brep_reader(brep_string, indx):

    shape_set = BRepTools_ShapeSet()
    shape_set.ReadFromString(brep_string)
    shape = shape_set.Shape(shape_set.NbShapes())
    location = shape_set.Locations().Location(indx)
    shape.Location(location)

    stl_writer = StlAPI_Writer()

    with siconos.io.tools.tmpfile(suffix=".stl") as tmpf:
        stl_writer.Write(shape, tmpf[1])
        tmpf[0].flush()

        reader = vtk.vtkSTLReader()
        reader.SetFileName(tmpf[1])
        reader.Update()

        return reader


def occ_load_brep(shape_ref, brep_class):
    """
    Load a BRep shape from a HDF5 dataset

    Parameters
    ----------
    shape_ref : h5py.Dataset
        HDF5 dataset containing either a raw BRep string or an encoded
        shape set with 'occ_indx'.

    brep_class : type
        Class used to instantiate the shape.
        Defaults to `siconos.mechanics.collision.occ.OccContactShape`.

    Returns
    -------
    brep : shape_class instance
        Loaded OpenCascade shape wrapped in OccContactShape)
    """
    if "occ_indx" in shape_ref.attrs:
        shape_set = BRepTools_ShapeSet()
        shape_set.ReadFromString(shape_ref[:][0])
        the_shape = shape_set.Shape(shape_set.NbShapes())
        location = shape_set.Locations().Location(shape_ref.attrs["occ_indx"])
        the_shape.Location(location)
        brep = brep_class()
        brep.setData(the_shape)

    else:
        # raw brep
        brep = brep_class()
        brep.importBRepFromString(shape_ref[:][0])

    return brep
