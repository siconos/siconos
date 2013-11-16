#
# from vtk files to bullet shapes
#

import os
import shlex
import vtk
from Siconos.Mechanics.ContactDetection.Bullet import btVector3, \
    btConvexHullShape, btCylinderShape, btBoxShape, btSphereShape, \
    btConeShape, btCapsuleShape, btCompoundShape, btTriangleIndexVertexArray, btGImpactMeshShape

import numpy as np

#
# load .vtp file
#
def loadShape(shape_filename):
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
        apoints = np.empty((num_points,3))
        for i in range(0, points.GetNumberOfTuples()):
            p = points.GetTuple(i)
            apoints[i,0] = p[0]
            apoints[i,1] = p[1]
            apoints[i,2] = p[2]

        aindices = np.empty((num_triangles,3), dtype=np.int32)

        for i in range(0, num_triangles):
            c = polydata.GetCell(i)
            aindices[i,0] = c.GetPointIds().GetId(0)
            aindices[i,1] = c.GetPointIds().GetId(1)
            aindices[i,2] = c.GetPointIds().GetId(2)

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

    def __init__(self, ref_filename='ref.txt'):
        self._ref_filename = ref_filename
        self._url = list()
        self._attributes = list()
        self._shape = dict()
        self._tri = dict()

        with open(self._ref_filename, 'r') as ref_file:
            for shape_url_line in ref_file:
                line_tokens = shlex.split(shape_url_line)
                shape_url = line_tokens[0]
                shape_attributes = [float(x) for x in (line_tokens[1:])]
                self._url.append(shape_url)
                self._attributes.append(shape_attributes)

        self._primitive = {'Cylinder': btCylinderShape,
                           'Sphere': btSphereShape,
                           'Box': btBoxShape,
                           'Cone': btConeShape,
                           'Compound': btCompoundShape,
                           'Capsule': btCapsuleShape}

    def at_index(self, index):

        if not index in self._shape:
            # load shape if it is an existing file
            if os.path.exists(self._url[index]):
                self._tri[index], self._shape[index] = loadShape(
                    self._url[index])
            else:
                # it must be a primitive with attributes
                name = self._url[index]
                primitive = self._primitive[name]
                attrs = self._attributes[index]
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
