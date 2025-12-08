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

class NativeSegmentShape(NativeShape):
    def __init__(self, x1, y1, x2, y2):
        self.params = [x1, y1, x2, y2]

class NativeBox2dShape(NativeShape):
    def __init__(self, a, b):
        self.params = [a, b]


def load_vtp_file(shape_ref):
    """
    Read a VTP shape from an HDF5 reference and return a loaded mesh.

    Parameters:
    ----------

    shape_ref: h5py.Dataset
        HDF5 dataset containing the VTP shape, stored as a string

    Returns
    -------
    mesh : SiconosMesh or SiconosConvexHull
        Loaded mesh object

    """
    assert shape_ref.attrs["type"] == "vtp"
    if shape_ref.dtype != h5py_vlen_dtype(str):
        raise AssertionError(
            "load_vtp_file(), wrong dtype, must be h5py_vlen_dtype(str) "
        )

    with siconos.io.tools.tmpfile() as tmpf:
        data = shape_ref[:][0]
        # fix compatibility with h5py version
        # to be removed in the future
        tmpf[0].write(data.decode("utf-8"))
        tmpf[0].flush()
        scale = shape_ref.attrs.get("scale", None)
        mesh, dims = siconos.io.tools.load_siconos_mesh(tmpf[1], scale=scale)
        mesh.setInsideMargin(shape_ref.attrs.get("insideMargin", min(dims) * 0.02))
        mesh.setOutsideMargin(shape_ref.attrs.get("outsideMargin", 0))
    return mesh


def load_heightmap(hm_data):
    # a heightmap defined by a matrix
    r = hm_data.attrs["rect"]
    if len(r) != 2:
        raise AssertionError("len(r) != 2")
    # assert(len(r) == 2)
    hm = siconos.mechanics.collision.SiconosHeightMap(np.array(hm_data, dtype=np.float64, order="F"), r[0], r[1])
    # dims = list(r) + [np.max(hm_data) - np.min(hm_data)]
    # hm.setInsideMargin(
    #    hm_data.attrs.get('insideMargin', np.min(dims) * 0.02))
    hm.setInsideMargin(hm_data.attrs.get("insideMargin", 0))
    hm.setOutsideMargin(hm_data.attrs.get("outsideMargin", 0))
    return hm


def load_convex(shape_ref, dimension):

    points = np.array(shape_ref, dtype=np.float64, order="F")
    if dimension == 3:
        convex = siconos.mechanics.collision.SiconosConvexHull(points)
        dims = [
            points[:, 0].max() - points[:, 0].min(),
            points[:, 1].max() - points[:, 1].min(),
            points[:, 2].max() - points[:, 2].min(),
        ]
    elif dimension == 2:
        convex = siconos.mechanics.collision.SiconosConvexHull2d(points)
        dims = [
            points[:, 0].max() - points[:, 0].min(),
            points[:, 1].max() - points[:, 1].min(),
        ]

    avoid_internal_edge_contact = shape_ref.attrs.get(
        "avoid_internal_edge_contact", False
    )
    if avoid_internal_edge_contact:
        convex.setAvoidInternalEdgeContact(True)

    convex.setInsideMargin(shape_ref.attrs.get("insideMargin", min(dims) * 0.02))
    convex.setOutsideMargin(shape_ref.attrs.get("outsideMargin", 0))
    return convex


class ShapeCollection:
    """
    Instantiation of added contact shapes

    Parameters
    ----------

    io: reference object. Either a MechanicsHdf5Runner
        or a file name (vtk file describing a mesh)
    collision_margin: tolerance for collision, defaults to 0.04
    backend: chosen backend, defaults to bullet

    """

    # Note FP: can't work with vtk input, can it?

    def __init__(self, io, collision_margin=None, backend="bullet"):
        self._io = io
        self._shapes = dict()
        self._tri = dict()
        if collision_margin is None:
            collision_margin = 0.04
        self._collision_margin = collision_margin
        self.backend = backend
        if backend == "vnative":
            self._primitive = {
                "Disk": NativeDiskShape,
                "Circle": NativeCircleShape,
                "Line": NativeLineShape,
                "Segment": NativeSegmentShape,
                "Box2d": NativeBox2dShape,
            }
        elif backend == "native":
            self._primitive = {
                "Disk": NativeDiskShape,
                "Circle": NativeCircleShape,
                "Line": NativeLineShape,
            }
        elif backend == "bullet":
            self._primitive = {
                "Sphere": siconos.mechanics.collision.SiconosSphere,
                "Box": siconos.mechanics.collision.SiconosBox,
                "Cylinder": siconos.mechanics.collision.SiconosCylinder,
                "Capsule": siconos.mechanics.collision.SiconosCapsule,
                "Cone": siconos.mechanics.collision.SiconosCone,
                "Plane": siconos.mechanics.collision.SiconosPlane,
                "Disk": siconos.mechanics.collision.SiconosDisk,
                "Segment": NativeSegmentShape, # transformed into Box2d
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

            shape_ref = self.shape(shape_name)

            # load shape if it is an existing file
            if (
                not isinstance(self.url(shape_name), str)
                and "primitive" not in shape_ref.attrs
            ):
                shape_type = shape_ref.attrs["type"]
                # --- case 1: load a vtp file ---
                # assume a vtp file (xml) stored in a string buffer
                if shape_type == "vtp":
                    self._shapes[shape_name] = load_vtp_file(shape_ref)

                # --- case 2: load a step file ---
                # occ only
                elif shape_type in ["step", "stp", "iges", "igs"]:
                    assert self.backend == "occ"
                    import siconos.io.occ_tools

                    if shape_ref.dtype != h5py_vlen_dtype(str):
                        raise AssertionError(
                            "ShapeCollection, get: dtype must be h5py_vlen_dtype(str)"
                        )

                    content = shape_ref[:][0].decode("utf-8")
                    with siconos.io.tools.tmpfile(
                        contents=content, suffix=".{0}".format(shape_type)
                    ) as tmpf:
                        comp = siconos.io.occ_tools.occ_load_file(tmpf[1])
                    self._shapes[shape_name] = comp
                    self._io._keep.append(self._shapes[shape_name])

                # --- case 3: load brep  ---
                elif shape_type in ["brep"]:
                    assert self.backend == "occ"
                    if shape_class is None:
                        brep_class = siconos.mechanics.collision.occ.OccContactShape
                    else:
                        brep_class = shape_class

                    if "contact" not in shape_ref.attrs:
                        self._shapes[shape_name] = siconos.io.occ_tools.occ_load_brep(
                            shape_ref, brep_class
                        )
                        self._io._keep.append(self._shapes[shape_name])

                    else:
                        required_keys = {"contact", "contact_index", "brep"}
                        missing = required_keys - shape_ref.attrs.keys()
                        if missing:
                            raise AssertionError(
                                f"Missing required attribute(s): {', '.join(missing)}"
                            )

                        contact_index = shape_ref.attrs["contact_index"]

                        ref_brep = self.get(shape_ref.attrs["brep"], shape_class)
                        if shape_ref.attrs["contact"] == "Face":
                            if face_class is None:
                                face_maker = (
                                    siconos.mechanics.collision.occ.OccContactFace
                                )
                            else:
                                face_maker = face_class

                            self._shapes[shape_name] = face_maker(
                                brep_class(ref_brep), contact_index
                            )

                        elif shape_ref.attrs["contact"] == "Edge":
                            if edge_class is None:
                                edge_maker = (
                                    siconos.mechanics.collision.occ.OccContactEdge
                                )
                            else:
                                edge_maker = edge_class
                            self._shapes[shape_name] = edge_maker(
                                ref_brep, contact_index
                            )

                        self._io._keep.append(self._shapes[shape_name])

                # --- case 5: load heightmap ---
                elif shape_type in ["heightmap"]:
                    self._shapes[shape_name] = load_heightmap(shape_ref)

                # --- case 6: load convex ---
                elif shape_type in ["convex"]:
                    # a convex point set
                    self._shapes[shape_name] = load_convex(
                        shape_ref, self._io._dimension
                    )

            elif isinstance(self.url(shape_name), str) and os.path.exists(
                self.url(shape_name)
            ):
                self._tri[shape_name], self._shapes[shape_name] = loadMesh(
                    self.url(shape_name), self._collision_margin
                )
                # Warning FP: can't find the function loadMesh anywhere?
                # Where does it supposed to come from?

            else:

                # it must be a primitive with attributes
                if isinstance(self.url(shape_name), str):
                    name = self.url(shape_name)
                    attrs = [float(x) for x in shape_ref[0]]
                else:
                    name = shape_ref.attrs["primitive"]
                    attrs = [float(x) for x in shape_ref[0]]
                primitive = self._primitive[name]

                if name in ["Box"]:
                    attrs = np.array(attrs, dtype=np.float64, order="F")
                    box = primitive(attrs)
                    self._shapes[shape_name] = box
                    box.setInsideMargin(
                        shape_ref.attrs.get("insideMargin", min(attrs) * 0.02)
                    )
                    box.setOutsideMargin(shape_ref.attrs.get("outsideMargin", 0))
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
                    shp = shape_ref
                    prim.setInsideMargin(
                        shp.attrs.get("insideMargin", min(attrs) * 0.02)
                    )
                    prim.setOutsideMargin(shp.attrs.get("outsideMargin", 0))

        return self._shapes[shape_name]
