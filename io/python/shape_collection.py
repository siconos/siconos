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
#

"""Class and tools to handle shapes in siconos.io.mechanics_run"""

import os
import numpy as np
import h5py
import siconos.mechanics.collision
import siconos.io.tools


# fix compatibility with h5py version
if hasattr(h5py, "vlen_dtype"):
    h5py_vlen_dtype = h5py.vlen_dtype
elif hasattr(h5py, "new_vlen"):
    h5py_vlen_dtype = h5py.new_vlen


class NativeShape:

    def setInsideMargin(self, m):
        self.insidemargin = m

    def setOutsideMargin(self, m):
        self.outsideMargin = m


class NativeDiskShape(NativeShape):
    def __init__(self, radius):
        self.radius = radius


class NativeCircleShape(NativeShape):
    def __init__(self, radius):
        self.radius = radius


class NativeLineShape(NativeShape):
    def __init__(self, a, b, c):
        self.params = [a, b, c]


class ShapeCollection:
    """
    Instantiation of added contact shapes
    """

    def __init__(self, io, collision_margin=0.04, backend="bullet"):
        self._io = io
        self._shapes = dict()
        self._tri = dict()

        if backend == "native":
            self._primitive = {
                "Disk": NativeDiskShape,
                "Circle": NativeCircleShape,
                "Line": NativeLineShape,
            }
        else:
            self._primitive = {
                "Sphere": siconos.mechanics.collision.SiconosSphere,
                "Box": siconos.mechanics.collision.SiconosBox,
                "Cylinder": siconos.mechanics.collision.SiconosCylinder,
                "Capsule": siconos.mechanics.collision.SiconosCapsule,
                "Cone": siconos.mechanics.collision.SiconosCone,
                "Plane": siconos.mechanics.collision.SiconosPlane,
                "Disk": siconos.mechanics.collision.SiconosDisk,
                "Box2d": siconos.mechanics.collision.SiconosBox2d,
            }

    def shape(self, shape_name):
        return self._io.shapes()[shape_name]

    def attributes(self, shape_name):
        return self._io.shapes()[shape_name].attrs

    def url(self, shape_name):
        if "url" in self.attributes(shape_name):
            shape_url = self.shape(shape_name).attrs["url"]

        elif "filename" in self.attributes(shape_name):
            shape_url = self.shape(shape_name).attrs["filename"]

        else:
            shape_url = self.shape(shape_name)

        return shape_url

    def get(
        self,
        shape_name,
        shape_class=None,
        face_class=None,
        edge_class=None,
        new_instance=False,
    ):

        if new_instance or shape_name not in self._shapes:

            # load shape if it is an existing file
            if not isinstance(
                self.url(shape_name), str
            ) and "primitive" not in self.attributes(shape_name):
                # assume a vtp file (xml) stored in a string buffer

                if self.attributes(shape_name)["type"] == "vtp":
                    if self.shape(shape_name).dtype == h5py_vlen_dtype(str):
                        with siconos.io.tools.tmpfile() as tmpf:
                            data = self.shape(shape_name)[:][0]
                            # fix compatibility with h5py version
                            # to be removed in the future
                            tmpf[0].write(data.decode("utf-8"))
                            tmpf[0].flush()
                            scale = self.attributes(shape_name).get("scale", None)
                            mesh, dims = siconos.io.tools.load_siconos_mesh(
                                tmpf[1], scale=scale
                            )
                            self._shapes[shape_name] = mesh
                            mesh.setInsideMargin(
                                self.shape(shape_name).attrs.get(
                                    "insideMargin", min(dims) * 0.02
                                )
                            )
                            mesh.setOutsideMargin(
                                self.shape(shape_name).attrs.get("outsideMargin", 0)
                            )
                    else:
                        raise AssertionError(
                            'ShapeCollection.get(), attributes(shape_name)["type"] != vtp '
                        )

                elif self.attributes(shape_name)["type"] in ["step", "stp"]:
                    from OCC.Core.STEPControl import STEPControl_Reader
                    from OCC.Core.BRep import BRep_Builder
                    from OCC.Core.TopoDS import TopoDS_Compound
                    from OCC.Core.IFSelect import (
                        IFSelect_RetDone,
                        IFSelect_ItemsByEntity,
                    )

                    builder = BRep_Builder()
                    comp = TopoDS_Compound()
                    builder.MakeCompound(comp)

                    if self.shape(shape_name).dtype != h5py_vlen_dtype(str):
                        raise AssertionError(
                            "self.shape(shape_name).dtype != h5py_vlen_dtype(str)"
                        )

                    tmp_contents = self.shape(shape_name)[:][0].decode("utf-8")

                    with siconos.io.tools.tmpfile(contents=tmp_contents) as tmpf:
                        step_reader = STEPControl_Reader()

                        status = step_reader.ReadFile(tmpf[1])

                        if status == IFSelect_RetDone:  # check status
                            failsonly = False
                            step_reader.PrintCheckLoad(
                                failsonly, IFSelect_ItemsByEntity
                            )
                            step_reader.PrintCheckTransfer(
                                failsonly, IFSelect_ItemsByEntity
                            )

                            #  ok = step_reader.TransferRoot(1)
                            step_reader.TransferRoots()
                            # VA : We decide to loads all shapes in the step file
                            nbs = step_reader.NbShapes()

                            for i in range(1, nbs + 1):
                                shape = step_reader.Shape(i)
                                builder.Add(comp, shape)

                            self._shapes[shape_name] = comp
                            self._io._keep.append(self._shapes[shape_name])

                elif self.attributes(shape_name)["type"] in ["iges", "igs"]:
                    from OCC.Core.IGESControl import IGESControl_Reader
                    from OCC.Core.BRep import BRep_Builder
                    from OCC.Core.TopoDS import TopoDS_Compound
                    from OCC.Core.IFSelect import (
                        IFSelect_RetDone,
                        IFSelect_ItemsByEntity,
                    )

                    builder = BRep_Builder()
                    comp = TopoDS_Compound()
                    builder.MakeCompound(comp)

                    if not (self.shape(shape_name).dtype == h5py_vlen_dtype(str)):
                        raise AssertionError("ShapeCollection.get()")

                    # assert(self.shape(shape_name).dtype == h5py_vlen_dtype(str))

                    with tmpfile(contents=self.shape(shape_name)[:][0]) as tmpf:
                        iges_reader = IGESControl_Reader()

                        status = iges_reader.ReadFile(tmpf[1])

                        if status == IFSelect_RetDone:  # check status
                            failsonly = False
                            iges_reader.PrintCheckLoad(
                                failsonly, IFSelect_ItemsByEntity
                            )
                            iges_reader.PrintCheckTransfer(
                                failsonly, IFSelect_ItemsByEntity
                            )

                            iges_reader.TransferRoots()
                            nbs = iges_reader.NbShapes()

                            for i in range(1, nbs + 1):
                                shape = iges_reader.Shape(i)
                                builder.Add(comp, shape)

                            self._shapes[shape_name] = comp
                            self._io._keep.append(self._shapes[shape_name])

                elif self.attributes(shape_name)["type"] in ["brep"]:
                    if "contact" not in self.attributes(shape_name):

                        # the reference brep
                        if shape_class is None:
                            brep_class = occ.OccContactShape
                        else:
                            brep_class = shape_class

                        if "occ_indx" in self.attributes(shape_name):

                            from OCC.Core.BRepTools import BRepTools_ShapeSet

                            shape_set = BRepTools_ShapeSet()
                            shape_set.ReadFromString(self.shape(shape_name)[:][0])
                            the_shape = shape_set.Shape(shape_set.NbShapes())
                            location = shape_set.Locations().Location(
                                self.attributes(shape_name)["occ_indx"]
                            )
                            the_shape.Location(location)
                            brep = brep_class()
                            brep.setData(the_shape)

                        else:
                            # raw brep
                            brep = brep_class()
                            brep.importBRepFromString(self.shape(shape_name)[:][0])

                        self._shapes[shape_name] = brep
                        self._io._keep.append(self._shapes[shape_name])

                    else:
                        # a contact on a brep
                        if not ("contact" in self.attributes(shape_name)):
                            raise AssertionError(
                                "contact not in self.attributes(shape_name)"
                            )
                        if not ("contact_index" in self.attributes(shape_name)):
                            raise AssertionError(
                                "contact_index not in self.attributes(shape_name)"
                            )
                        if not ("brep" in self.attributes(shape_name)):
                            raise AssertionError(
                                "brep not in self.attributes(shape_name)"
                            )

                        # assert 'contact' in self.attributes(shape_name)
                        # assert 'contact_index' in self.attributes(shape_name)
                        # assert 'brep' in self.attributes(shape_name)

                        contact_index = self.attributes(shape_name)["contact_index"]

                        if shape_class is None:
                            brep_class = occ.OccContactShape
                        else:
                            brep_class = shape_class

                        ref_brep = self.get(
                            self.attributes(shape_name)["brep"], shape_class
                        )

                        if self.attributes(shape_name)["contact"] == "Face":
                            if face_class is None:
                                face_maker = occ.OccContactFace
                            else:
                                face_maker = face_class

                            self._shapes[shape_name] = face_maker(
                                brep_class(ref_brep), contact_index
                            )

                        elif self.attributes(shape_name)["contact"] == "Edge":
                            if edge_class is None:
                                edge_maker = occ.OccContactEdge
                            else:
                                edge_maker = edge_class
                            self._shapes[shape_name] = edge_maker(
                                ref_brep, contact_index
                            )

                        self._io._keep.append(self._shapes[shape_name])

                elif self.attributes(shape_name)["type"] in ["heightmap"]:
                    # a heightmap defined by a matrix
                    hm_data = self.shape(shape_name)
                    r = hm_data.attrs["rect"]
                    if len(r) != 2:
                        raise AssertionError("len(r) != 2")
                    # assert(len(r) == 2)
                    hm = siconos.mechanics.collision.SiconosHeightMap(
                        hm_data, r[0], r[1]
                    )
                    # dims = list(r) + [np.max(hm_data) - np.min(hm_data)]
                    # hm.setInsideMargin(
                    #    hm_data.attrs.get('insideMargin', np.min(dims) * 0.02))
                    hm.setInsideMargin(hm_data.attrs.get("insideMargin", 0))
                    hm.setOutsideMargin(hm_data.attrs.get("outsideMargin", 0))

                    self._shapes[shape_name] = hm

                elif self.attributes(shape_name)["type"] in ["convex"]:
                    # a convex point set
                    points = self.shape(shape_name)
                    points = np.array(points, dtype=np.float64, order="F")
                    if self._io._dimension == 3:
                        convex = siconos.mechanics.collision.SiconosConvexHull(points)
                        dims = [
                            points[:, 0].max() - points[:, 0].min(),
                            points[:, 1].max() - points[:, 1].min(),
                            points[:, 2].max() - points[:, 2].min(),
                        ]
                    elif self._io._dimension == 2:
                        convex = siconos.mechanics.collision.SiconosConvexHull2d(points)
                        dims = [
                            points[:, 0].max() - points[:, 0].min(),
                            points[:, 1].max() - points[:, 1].min(),
                        ]

                    avoid_internal_edge_contact = self.shape(shape_name).attrs.get(
                        "avoid_internal_edge_contact", False
                    )
                    if avoid_internal_edge_contact:
                        convex.setAvoidInternalEdgeContact(True)

                    convex.setInsideMargin(
                        self.shape(shape_name).attrs.get(
                            "insideMargin", min(dims) * 0.02
                        )
                    )
                    convex.setOutsideMargin(
                        self.shape(shape_name).attrs.get("outsideMargin", 0)
                    )
                    self._shapes[shape_name] = convex

            elif isinstance(self.url(shape_name), str) and os.path.exists(
                self.url(shape_name)
            ):
                self._tri[shape_name], self._shapes[shape_name] = loadMesh(
                    self.url(shape_name), _collision_margin
                )
            else:
                # it must be a primitive with attributes
                if isinstance(self.url(shape_name), str):
                    name = self.url(shape_name)
                    attrs = [float(x) for x in self.shape(shape_name)[0]]
                else:
                    name = self.attributes(shape_name)["primitive"]
                    attrs = [float(x) for x in self.shape(shape_name)[0]]
                primitive = self._primitive[name]

                if name in ["Box"]:
                    attrs = np.array(attrs, dtype=np.float64, order="F")
                    box = primitive(attrs)
                    self._shapes[shape_name] = box
                    box.setInsideMargin(
                        self.shape(shape_name).attrs.get(
                            "insideMargin", min(attrs) * 0.02
                        )
                    )
                    box.setOutsideMargin(
                        self.shape(shape_name).attrs.get("outsideMargin", 0)
                    )
                # elif name in ['Compound']:
                #     obj1 = attrs[0]
                #     orig1 = attrs[1:4]
                #     orie1 = attrs[4:8]
                #     obj2 = attrs[8]
                #     orig2 = attrs[9:12]
                #     orie2 = attrs[12:16]
                #     bcols = btCompoundShape()
                #     bcols.addChildShape(...
                else:  # e.g. name in ['Sphere']:
                    prim = self._shapes[shape_name] = primitive(*attrs)
                    shp = self.shape(shape_name)
                    prim.setInsideMargin(
                        shp.attrs.get("insideMargin", min(attrs) * 0.02)
                    )
                    prim.setOutsideMargin(shp.attrs.get("outsideMargin", 0))

        return self._shapes[shape_name]
