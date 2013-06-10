#
# from vtk files to bullet shapes
#

import os
import shlex
import vtk
from Siconos.Mechanics.ContactDetection.Bullet import btVector3, \
    btConvexHullShape


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


class Collection(dict):

    """
    Returns Bullet convex hull shapes from .vtp filenames given in a
    reference file

    The access is done like with dictionary with shape name (filename
    without extension) as key
    """

    def __init__(self, ref_filename = 'ref.txt'):
        self._ref_filename = ref_filename
        self._url = dict()
        self._shape = dict()
        self._ref = dict()
        with open(self._ref_filename, 'r') as ref_file:
            index = 0
            for shape_url_line in ref_file:
                shape_url = shlex.split(shape_url_line)[0]
                name = os.path.splitext(os.path.basename(shape_url))[0]
                self._url[name] = shape_url
                self._shape[index] = name
                self._ref[name] = index
                index += 1

        super(Collection, self).__init__()

    def __getitem__(self, name):

        if not name in self:
            # load shape
            self[name] = loadConvexHullShape(
                self._url[name])

        return super(Collection, self).__getitem__(name)

    def at_index(self, ind):
        return self[self._shape[ind]]

    def index(self, name):
        return self._ref[name]
