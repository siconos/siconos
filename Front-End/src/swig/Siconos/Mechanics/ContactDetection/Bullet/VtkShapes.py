#
# from vtk files to bullet shapes
#

import os
import shlex
import vtk
from Siconos.Mechanics.ContactDetection.Bullet import btVector3, \
    btConvexHullShape, btCylinderShape, btBoxShape, btSphereShape, \
    btConeShape, btCapsuleShape, btCompoundShape, btTriangleIndexVertexArray, \
    btGImpactMeshShape

import numpy as np

from operator import itemgetter
import tempfile
from contextlib import contextmanager

import h5py

#
# a context manager for a *named* temp file
#
@contextmanager
def tmpfile():
    (_, filename) = tempfile.mkstemp()
    fid = open(filename, 'w')
    yield (fid, filename)
    fid.close()
    os.remove(filename)

#
# load .vtp file
#
def loadMesh(shape_filename):
    """
    loads a vtk .vtp file and returns a Bullet concave shape
    WARNING triangles cells assumed!
    """

    reader = vtk.vtkXMLPolyDataReader()
    reader.SetFileName(shape_filename)
    reader.Update()

    polydata = reader.GetOutput()
    points = polydata.GetPoints().GetData()
    num_points = points.GetNumberOfTuples()
    num_triangles = polydata.GetNumberOfCells()

    print num_points, num_triangles

    keep = None
    shape = None

    if polydata.GetCellType(0) == 5:
        apoints = np.empty((num_points, 3))
        for i in range(0, points.GetNumberOfTuples()):
            p = points.GetTuple(i)
            apoints[i, 0] = p[0]
            apoints[i, 1] = p[1]
            apoints[i, 2] = p[2]

        aindices = np.empty((num_triangles, 3), dtype=np.int32)

        for i in range(0, num_triangles):
            c = polydata.GetCell(i)
            aindices[i, 0] = c.GetPointIds().GetId(0)
            aindices[i, 1] = c.GetPointIds().GetId(1)
            aindices[i, 2] = c.GetPointIds().GetId(2)

        tri = btTriangleIndexVertexArray(apoints, aindices)

        shape = btGImpactMeshShape(tri)
        shape.updateBound()

        keep = tri, apoints, aindices

    else:  # assume convex shape
        coors = dict()
        for i in range(0, points.GetNumberOfTuples()):
            coors[points.GetTuple(i)] = 1

        shape = btConvexHullShape()
        for p in coors:
                shape.addPoint(btVector3(*p))

    return keep, shape


class Collection():

    """
    collect Bullet primitives or convex hull shapes from .vtp
    filenames given in a reference file
    """

    def __init__(self, ref = 'ref.txt'):
        self._ref = ref
        self._url = list()
        self._attributes = list()
        self._shape = dict()
        self._tri = dict()

        if isinstance(self._ref, str):
            with open(self._ref, 'r') as ref_file:
                for shape_url_line in ref_file:
                    line_tokens = shlex.split(shape_url_line)
                    shape_url = line_tokens[0]
                    shape_attributes = [float(x) for x in (line_tokens[1:])]
                    self._url.append(shape_url)
                    self._attributes.append(shape_attributes)

        # assume hdf5 file
        else:
            acc = []
            for shape_name in self._ref['data']['ref']:
                shape_attributes = self._ref['data']['ref'][shape_name][:]
                shape_id = self._ref['data']['ref'][shape_name].attrs['id']
                if 'url' in self._ref['data']['ref'][shape_name].attrs:
                    shape_url = self._ref['data']['ref'][shape_name].\
                        attrs['url']

                elif 'filename' in self._ref['data']['ref'][shape_name].attrs:
                    shape_url = self._ref['data']['ref'][shape_name].\
                        attrs['filename']

                else:
                    shape_url = self._ref['data']['ref'][shape_name]

                acc.append((shape_attributes, shape_url, shape_id))

            sacc = sorted(acc, key=itemgetter(2))
            self._url = [u for a, u, i in sacc]
            self._attributes = [a for a, u, i in sacc]

        self._primitive = {'Cylinder': btCylinderShape,
                           'Sphere': btSphereShape,
                           'Box': btBoxShape,
                           'Cone': btConeShape,
                           'Compound': btCompoundShape,
                           'Capsule': btCapsuleShape}

    def at_index(self, index):

        if not index in self._shape:

            # load shape if it is an existing file
            if not isinstance(self._url[index], str) and \
                not 'primitive' in self._url[index].attrs:
                # assume a vtp file (xml) stored in a string buffer
                if self._url[index].dtype == h5py.new_vlen(str):
                    with tmpfile() as tmpf:
                        data = self._url[index][:][0]
                        tmpf[0].write(data)
                        tmpf[0].flush()
                        self._tri[index], self._shape[index] = loadMesh(
                            tmpf[1])
                else:
                    # a convex point set
                    convex = btConvexHullShape()
                    for points in self._url[index]:
                        convex.addPoint(btVector3(float(points[0]),
                                                  float(points[1]),
                                                  float(points[2])))
                    self._shape[index] = convex

            elif isinstance(self._url[index], str) and \
                os.path.exists(self._url[index]):
                self._tri[index], self._shape[index] = loadMesh(
                    self._url[index])
            else:
                # it must be a primitive with attributes
                if isinstance(self._url[index], str):
                    name = self._url[index]
                    attrs = [float(x) for x in self._attributes[index]]
                else:
                    name = self._url[index].attrs['primitive']
                    attrs = [float(x) for x in self._attributes[index][0]]
                primitive = self._primitive[name]

                if name in ['Box']:
                    self._shape[index] = primitive(btVector3(attrs[0] / 2,
                                                             attrs[1] / 2,
                                                             attrs[2] / 2))
                elif name in ['Cylinder']:
                    self._shape[index] = primitive(btVector3(attrs[0],
                                                             attrs[1] / 2,
                                                             attrs[0]))
                # elif name in ['Compound']:
                #     obj1 = attrs[0]
                #     orig1 = attrs[1:4]
                #     orie1 = attrs[4:8]
                #     obj2 = attrs[8]
                #     orig2 = attrs[9:12]
                #     orie2 = attrs[12:16]
                #     bcols = btCompoundShape()
                #     bcols.addChildShape(...
                else:
                    self._shape[index] = primitive(*attrs)

        return self._shape[index]
