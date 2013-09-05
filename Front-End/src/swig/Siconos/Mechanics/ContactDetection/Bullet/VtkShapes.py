#
# from vtk files to bullet shapes
#

import os
import shlex
import vtk
from Siconos.Mechanics.ContactDetection.Bullet import btVector3, \
    btConvexHullShape, btCylinderShape, btBoxShape, btSphereShape, \
    btConeShape, btCapsuleShape, btCompoundShape


#
# load .vtp file
#
def loadConvexHullShape(shape_filename):
    """
    loads a vtk .vtp file and returns a Bullet convex hull shape
    """
    reader = vtk.vtkXMLPolyDataReader()
    reader.SetFileName(shape_filename)
    reader.Update()
    polydata = reader.GetOutput()
    points = polydata.GetPoints().GetData()
    coors = dict()
    for i in range(0, points.GetNumberOfTuples()):
        coors[points.GetTuple(i)] = 1

    convex_hull_shape = btConvexHullShape()
    for p in coors:
        convex_hull_shape.addPoint(btVector3(*p))

    return convex_hull_shape


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
                self._shape[index] = loadConvexHullShape(
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
