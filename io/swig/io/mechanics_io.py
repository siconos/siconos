# Mechanics IO

from __future__ import print_function

import os
import sys

from math import cos, sin, asin, atan2
from scipy import constants

import numpy as np
import h5py
import bisect
import time
import pickle

import tempfile
from contextlib import contextmanager

# Siconos imports
import siconos.numerics as Numerics
from siconos.kernel import \
    EqualityConditionNSL, \
    Interaction, DynamicalSystem, TimeStepping
import siconos.kernel as Kernel

# Siconos Mechanics imports
from siconos.mechanics.collision.tools import Contactor, Volume
from siconos.mechanics import joints
from siconos.io.io_base import MechanicsIO
from siconos.io.FrictionContactTrace import FrictionContactTrace

# Currently we must choose between two implementations
use_original = False
use_proposed = True


def set_implementation(i):
    global use_original, use_proposed
    if i == 'original':
        use_original = have_original
        use_proposed = not have_original
        setup_default_classes()
        return use_original
    elif i == 'proposed':
        use_proposed = have_proposed
        use_original = not have_original
        setup_default_classes()
        return use_proposed
    return False

# For 'proposed' implementation, it is necessary to select a back-end,
# although currently only Bullet is supported for general objects.
backend = 'bullet'


def set_backend(b):
    global backend
    backend = b
    setup_default_classes()

have_proposed = False
have_original = False
have_bullet = False
have_occ = False

# Imports for 'proposed' implementation
try:
    from siconos.mechanics.collision import BodyDS, \
        SiconosSphere, SiconosBox, SiconosCylinder, SiconosPlane, \
        SiconosConvexHull, SiconosContactor, SiconosContactorSet, \
        SiconosMesh, SiconosHeightMap

    try:
        from siconos.mechanics.collision.bullet import \
            SiconosBulletCollisionManager, SiconosBulletOptions
        have_bullet = True
    except:
        have_bullet = False

    have_proposed = True
except:
    have_proposed = False
    use_proposed = False

# Imports for 'original' implementation
try:
    from siconos.mechanics.collision.bullet import \
        BulletDS, BulletWeightedShape, btScalarSize, \
        btCollisionObject, BulletTimeStepping, BulletSpaceFilter

    from siconos.mechanics.collision.bullet import btVector3, \
        btConvexHullShape, btCylinderShape, btBoxShape, btSphereShape, \
        btConeShape, btCapsuleShape, btCompoundShape, btTriangleIndexVertexArray, \
        btGImpactMeshShape
    have_bullet = True
    have_original = True

except:
    have_original = False
    use_original = False

# Shared Bullet imports
try:
    from siconos.mechanics.collision.bullet import \
        btScalarSize, btQuaternion, btTransform, \
        btVector3, quatRotate
    from siconos.mechanics.collision.bullet import \
        __mul__ as mul
except:
    have_bullet = False

# OCC imports
try:
    from siconos.mechanics import occ
    have_occ = True
except:
    have_occ = False

### Configuration

default_manager_class = None
default_simulation_class = None
default_body_class = None
use_bullet = False

def setup_default_classes():
    global default_manager_class
    global default_simulation_class
    global default_body_class
    global use_bullet
    if use_proposed:
        if backend == 'bullet':
            def m(model, options):
                if options is None:
                    options = SiconosBulletOptions()
                return SiconosBulletCollisionManager(options)
            default_manager_class = m
            use_bullet = have_bullet
        default_simulation_class = TimeStepping
        default_body_class = BodyDS
    elif use_original:
        if backend == 'bullet':
            default_manager_class = lambda model,options: BulletSpaceFilter(model)
            default_simulation_class = BulletTimeStepping
            default_body_class = BulletDS
            use_bullet = have_bullet
        elif backend == 'occ':
            default_manager_class = lambda model,options: occ.OccSpaceFilter(model)
            default_simulation_class = occ.OccTimeStepping
            default_body_class = occ.OccBody
            use_bullet = have_bullet

setup_default_classes()

### Utility functions

def floatv(v):
    return [float(x) for x in v]


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
def tmpfile(suffix='', prefix='siconos_io', contents=None):
    """
    A context manager for a named temporary file.
    """
    (_, tfilename) = tempfile.mkstemp(suffix=suffix, prefix=prefix)
    fid = open(tfilename, 'w')
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
    os.remove(tfilename)


class Timer():

    def __init__(self):
        self._t0 = time.clock()

    def elapsed(self):
        return time.clock() - self._t0

    def update(self):
        self._t0 = time.clock()


def warn(msg):
    sys.stderr.write('{0}: {1}'.format(sys.argv[0], msg))


def log(fun, with_timer=False):
    if with_timer:
        t = Timer()

        def logged(*args):
            t.update()
            print('{0} ...'.format(fun.__name__), end='')
            fun(*args)
            print('..... {0} s'.format(t.elapsed()))
        return logged
    else:
        def silent(*args):
            fun(*args)
        return silent


def object_id(obj):
    """returns an unique object identifier"""
    return obj.__hash__()


def apply_gravity(body):
    g = constants.g
    weight = [0, 0, - body.scalarMass() * g]
    body.setFExtPtr(weight)


def group(h, name, must_exist=True):
    try:
        return h[name]
    except KeyError:
        if must_exist:
            return h.create_group(name)
        else:
            try:
                return h.create_group(name)
            except ValueError:
                # could not create group, return None
                # (file is probably in read-only mode)
                return None

def data(h, name, nbcolumns, use_compression=False):
    try:
        return h[name]
    except KeyError:
        comp = use_compression and nbcolumns > 0
        return h.create_dataset(name, (0, nbcolumns),
                                maxshape=(None, nbcolumns),
                                chunks=[None,(4000,nbcolumns)][comp],
                                compression=[None,'gzip'][comp],
                                compression_opts=[None,9][comp])


def add_line(dataset, line):
    dataset.resize(dataset.shape[0] + 1, 0)
    dataset[dataset.shape[0] - 1, :] = line


def str_of_file(filename):
    with open(filename, 'r') as f:
        return str(f.read())


class Quaternion():

    def __init__(self, *args):
        import vtk
        self._vtkmath = vtk.vtkMath()
        self._data = vtk.vtkQuaternion[float](*args)

    def __mul__(self, q):
        r = Quaternion()
        self._vtkmath.MultiplyQuaternion(self._data, q._data, r._data)
        return r

    def __getitem__(self, i):
        return self._data[i]

    def conjugate(self):
        r = Quaternion((self[0], self[1], self[2], self[3]))
        r._data.Conjugate()
        return r

    def rotate(self, v):
        pv = Quaternion((0, v[0], v[1], v[2]))
        rv = self * pv * self.conjugate()
        # assert(rv[0] == 0)
        return [rv[1], rv[2], rv[3]]

    def axisAngle(self):
        r = [0, 0, 0]
        a = self._data.GetRotationAngleAndAxis(r)
        return r, a


def phi(q0, q1, q2, q3):
    """
    Euler angle phi from quaternion.
    """
    return atan2(2*(q0*q1+q2*q3), 1-2*(q1*q1+q2*q2))


def theta(q0, q1, q2, q3):
    """
    Euler angle theta from quaternion.
    """
    return asin(2*(q0*q2-q3*q1))


def psi(q0, q1, q2, q3):
    """
    Euler angle psi from quaternion.
    """
    return atan2(2*(q0*q3+q1*q2), 1-2*(q2*q2+q3*q3))

# vectorized versions
phiv = np.vectorize(phi)
thetav = np.vectorize(theta)
psiv = np.vectorize(psi)


#
# load .vtp file
#
def loadMesh(shape_filename, collision_margin, scale=None):
    """
    loads a vtk .vtp file and returns a Bullet concave shape
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

    keep = None
    shape = None

    if polydata.GetCellType(0) == 5:
        apoints = np.empty((num_points, 3), dtype={4:'f4',8:'f8'}[btScalarSize()])
        for i in range(0, points.GetNumberOfTuples()):
            p = points.GetTuple(i)
            apoints[i, 0] = p[0]
            apoints[i, 1] = p[1]
            apoints[i, 2] = p[2]

        if scale is not None:
            apoints *= scale

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
        shape.setMargin(collision_margin)
        for p in coors:
                shape.addPoint(btVector3(*p))

    return keep, shape

def loadSiconosMesh(shape_filename, scale=None):
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

    keep = None
    shape = None

    if polydata.GetCellType(0) == 5:
        apoints = np.empty((num_points, 3), dtype={4:'f4',8:'f8'}[btScalarSize()])
        for i in range(0, points.GetNumberOfTuples()):
            p = points.GetTuple(i)
            apoints[i, 0] = p[0]
            apoints[i, 1] = p[1]
            apoints[i, 2] = p[2]

        if scale is not None:
            apoints *= scale

        aindices = np.empty(num_triangles*3, dtype=np.int)

        for i in range(0, num_triangles):
            c = polydata.GetCell(i)
            aindices[i*3 + 0] = c.GetPointIds().GetId(0)
            aindices[i*3 + 1] = c.GetPointIds().GetId(1)
            aindices[i*3 + 2] = c.GetPointIds().GetId(2)

        shape = SiconosMesh(list(aindices), apoints)
        dims = apoints.max(axis=0) - apoints.min(axis=0)

    else:  # assume convex shape
        coors = dict()
        for i in range(0, points.GetNumberOfTuples()):
            coors[points.GetTuple(i)] = 1
        coors = np.array(coors.keys())
        dims = coors.max(axis=0) - coors.min(axis=0)
        shape = SiconosConvexHull(coors.keys())

    return shape, dims

class ShapeCollection():

    """
    Instantiation of added contact shapes
    """

    def __init__(self, io, collision_margin=0.04):
        self._io = io
        self._shapes = dict()
        self._tri = dict()
        self._collision_margin=collision_margin
        # print('self._collision_margin',self._collision_margin)
        if use_proposed:

            self._primitive = {'Sphere': SiconosSphere,
                               'Box': SiconosBox,
                               'Cylinder': SiconosCylinder,
                               'Plane': SiconosPlane}

        elif use_original and use_bullet:

            self._primitive = {'Cylinder': btCylinderShape,
                               'Sphere': btSphereShape,
                               'Box': btBoxShape,
                               'Cone': btConeShape,
                               'Compound': btCompoundShape,
                               'Capsule': btCapsuleShape}
        else:
            self._primitive = dict()

    def shape(self, shape_name):
        return self._io.shapes()[shape_name]

    def attributes(self, shape_name):
        return self._io.shapes()[shape_name].attrs

    def url(self, shape_name):
        if 'url' in self.attributes(shape_name):
            shape_url = self.shape(shape_name).\
                attrs['url']

        elif 'filename' in self.attributes(shape_name):
            shape_url = self.shape(shape_name).\
                attrs['filename']

        else:
            shape_url = self.shape(shape_name)

        return shape_url

    def get(self, shape_name, shape_class=None, face_class=None,
            edge_class=None, new_instance=False):

        if new_instance or not shape_name in self._shapes:

            # load shape if it is an existing file
            if not isinstance(self.url(shape_name), str) and \
               not 'primitive' in self.attributes(shape_name):
                # assume a vtp file (xml) stored in a string buffer

                if self.attributes(shape_name)['type'] == 'vtp':
                    if self.shape(shape_name).dtype == h5py.new_vlen(str):
                        with tmpfile() as tmpf:
                            data = self.shape(shape_name)[:][0]
                            tmpf[0].write(data)
                            tmpf[0].flush()
                            scale = None
                            if 'scale' in self.attributes(shape_name):
                                scale = self.attributes(shape_name)['scale']
                            if use_proposed:
                                mesh, dims = loadSiconosMesh(tmpf[1], scale=scale)
                                self._shapes[shape_name] = mesh
                                mesh.setInsideMargin(
                                    self.shape(shape_name).attrs.get('insideMargin',
                                                                     min(dims)*0.02))
                                mesh.setOutsideMargin(
                                    self.shape(shape_name).attrs.get('outsideMargin',0))
                            elif use_original:
                                (self._tri[shape_name],
                                 self._shapes[shape_name]) = loadMesh(
                                     tmpf[1], self._collision_margin, scale=scale)
                    else:
                        assert False
                elif self.attributes(shape_name)['type'] in['step', 'stp']:
                    from OCC.STEPControl import STEPControl_Reader
                    from OCC.BRep import BRep_Builder
                    from OCC.TopoDS import TopoDS_Compound
                    from OCC.IFSelect import IFSelect_RetDone,\
                        IFSelect_ItemsByEntity

                    builder = BRep_Builder()
                    comp = TopoDS_Compound()
                    builder.MakeCompound(comp)

                    assert self.shape(shape_name).dtype == h5py.new_vlen(str)

                    with tmpfile(contents=self.shape(shape_name)[:][0]) as tmpf:
                        step_reader = STEPControl_Reader()

                        status = step_reader.ReadFile(tmpf[1])

                        if status == IFSelect_RetDone:  # check status
                            failsonly = False
                            step_reader.PrintCheckLoad(
                                failsonly, IFSelect_ItemsByEntity)
                            step_reader.PrintCheckTransfer(
                                failsonly, IFSelect_ItemsByEntity)

                            ok = step_reader.TransferRoot(1)
                            nbs = step_reader.NbShapes()

                            for i in range(1, nbs + 1):
                                shape = step_reader.Shape(i)
                                builder.Add(comp, shape)

                            self._shapes[shape_name] = comp
                            self._io._keep.append(self._shapes[shape_name])

                elif self.attributes(shape_name)['type'] in['brep']:
                    if not 'contact' in self.attributes(shape_name):

                        # the reference brep
                        if shape_class is None:
                            brep_class = occ.OccContactShape
                        else:
                            brep_class = shape_class

                        if 'occ_indx' in self.attributes(shape_name):

                            from OCC.BRepTools import BRepTools_ShapeSet
                            shape_set = BRepTools_ShapeSet()
                            shape_set.ReadFromString(
                                self.shape(shape_name)[:][0])
                            the_shape = shape_set.Shape(shape_set.NbShapes())
                            location = shape_set.Locations().Location(
                                self.attributes(shape_name)['occ_indx'])
                            the_shape.Location(location)
                            brep = brep_class()
                            brep.setData(the_shape)

                        else:
                            # raw brep
                            brep = brep_class()
                            brep.importBRepFromString(
                                self.shape(shape_name)[:][0])

                        self._shapes[shape_name] = brep
                        self._io._keep.append(self._shapes[shape_name])

                    else:
                        # a contact on a brep
                        assert 'contact' in self.attributes(shape_name)
                        assert 'contact_index' in self.attributes(shape_name)
                        assert 'brep' in self.attributes(shape_name)
                        contact_index = self.attributes(shape_name)['contact_index']

                        if shape_class is None:
                            brep_class = occ.OccContactShape
                        else:
                            brep_class = shape_class

                        ref_brep = self.get(
                            self.attributes(shape_name)['brep'], shape_class)

                        if self.attributes(shape_name)['contact'] == 'Face':
                            if face_class is None:
                                face_maker = occ.OccContactFace
                            else:
                                face_maker = face_class

                            self._shapes[shape_name] = \
                                face_maker(brep_class(ref_brep),
                                           contact_index)

                        elif self.attributes(shape_name)['contact'] == 'Edge':
                            if edge_class is None:
                                edge_maker = occ.OccContactEdge
                            else:
                                edge_maker = edge_class
                            self._shapes[shape_name] = \
                                edge_maker(ref_brep,
                                           contact_index)

                        self._io._keep.append(self._shapes[shape_name])

                elif self.attributes(shape_name)['type'] in ['heightmap']:

                    if use_proposed:
                        hm_data = self.shape(shape_name)
                        r = hm_data.attrs['rect']
                        assert(len(r)==2)
                        hm = SiconosHeightMap(hm_data, r[0], r[1])
                        dims = list(r) + [np.max(hm_data)-np.min(hm_data)]
                        hm.setInsideMargin(
                            hm_data.attrs.get('insideMargin', np.min(dims)*0.02))
                        hm.setOutsideMargin(
                            hm_data.attrs.get('outsideMargin', 0))

                        self._shapes[shape_name] = hm
                    else:
                        throw

                elif self.attributes(shape_name)['type'] in ['convex']:
                    # a convex point set
                    if use_proposed:
                        points = self.shape(shape_name)
                        convex = SiconosConvexHull(points)
                        dims = [points[:,0].max() - points[:,0].min(),
                                points[:,1].max() - points[:,1].min(),
                                points[:,2].max() - points[:,2].min()]
                        convex.setInsideMargin(
                            self.shape(shape_name).attrs.get('insideMargin',
                                                             min(dims)*0.02))
                        convex.setOutsideMargin(
                            self.shape(shape_name).attrs.get('outsideMargin', 0))
                    elif use_original and use_bullet:
                        convex = btConvexHullShape()
                        convex.setMargin(self._collision_margin)
                        for points in self.shape(shape_name):
                            convex.addPoint(btVector3(float(points[0]),
                                                      float(points[1]),
                                                      float(points[2])))
                    else:
                        throw
                    self._shapes[shape_name] = convex

                else:
                    throw

            elif isinstance(self.url(shape_name), str) and \
                    os.path.exists(self.url(shape_name)):
                self._tri[shape_name], self._shapes[shape_name] = loadMesh(
                    self.url(shape_name), _collision_margin)
            else:
                # it must be a primitive with attributes
                if isinstance(self.url(shape_name), str):
                    name = self.url(shape_name)
                    attrs = [float(x) for x in self.shape(shape_name)[0]]
                else:
                    name = self.attributes(shape_name)['primitive']
                    attrs = [float(x) for x in self.shape(shape_name)[0]]
                primitive = self._primitive[name]

                if name in ['Box']:
                    if use_proposed:
                        box = primitive(attrs)
                        self._shapes[shape_name] = box
                        box.setInsideMargin(
                            self.shape(shape_name).attrs.get('insideMargin',
                                                             min(attrs)*0.02))
                        box.setOutsideMargin(
                            self.shape(shape_name).attrs.get('outsideMargin', 0))
                    elif use_original and use_bullet:
                        self._shapes[shape_name] = primitive(
                            btVector3(attrs[0] / 2,
                                      attrs[1] / 2,
                                      attrs[2] / 2))

                elif name in ['Cylinder'] and not use_proposed:
                    self._shapes[shape_name] = primitive(btVector3(attrs[0],
                                                                   attrs[1]/2,
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
                else: # e.g. name in ['Sphere']:
                    prim = self._shapes[shape_name] = primitive(*attrs)
                    shp = self.shape(shape_name)
                    if use_proposed:
                        prim.setInsideMargin(
                            shp.attrs.get('insideMargin', min(attrs)*0.02))
                        prim.setOutsideMargin(shp.attrs.get('outsideMargin', 0))

        return self._shapes[shape_name]


class Hdf5():

    """a Hdf5 context manager reads at instantiation the translations and
       orientations of collision objects from hdf5 file

       It provides functions to output translations and orientations in
       the same file during simulation (output is done by default in
       pos.dat)

       with:
         time : float
         object_id : the object id (int)
         px, py, pz : components of the translation (float)
         ow, ox, oy oz : components of an unit quaternion (float)

    """

    def __init__(self, io_filename=None, mode='w',
                 broadphase=None, model=None, osi=None, shape_filename=None,
                 set_external_forces=None, gravity_scale=None, collision_margin=None,
                 use_compression=False, output_domains=False):

        if io_filename is None:
            self._io_filename = '{0}.hdf5'.format(
                os.path.splitext(os.path.basename(sys.argv[0]))[0])
        else:
            self._io_filename = io_filename
        self._mode = mode
        self._broadphase = broadphase
        self._model = model
        self._osi = osi
        self._static = {}
        self._shape = None
        self._shapeid = dict()
        self._pinterid = dict()
        self._static_data = None
        self._velocities_data = None
        self._dynamic_data = None
        self._cf_data = None
        self._domain_data = None
        self._solv_data = None
        self._input = None
        self._nslaws_data = None
        self._nslaws = dict()
        self._out = None
        self._data = None
        self._ref = None
        self._permanent_interactions = None
        self._occ_contactors = dict()
        self._joints = None
        self._boundary_conditions = None
        self._io = MechanicsIO()
        self._set_external_forces = set_external_forces
        self._shape_filename = shape_filename
        self._number_of_shapes = 0
        self._number_of_permanent_interactions = 0
        self._number_of_dynamic_objects = 0
        self._number_of_static_objects = 0
        self._gravity_scale = gravity_scale
        self._collision_margin = collision_margin
        self._output_frequency = 1
        self._keep = []
        self._scheduled_births = []
        self._births = dict()
        self._initializing = True
        self._use_compression = use_compression
        self._should_output_domains = output_domains
        self._contact_index_set = 1

    def __enter__(self):
        if self._set_external_forces is None:
            self._set_external_forces = self.apply_gravity

        if self._gravity_scale is None:
            self._gravity_scale = 1  # 1 => m, 1/100. => cm

        self._out = h5py.File(self._io_filename, self._mode)
        self._data = group(self._out, 'data')
        self._ref = group(self._data, 'ref')
        self._permanent_interactions = group(self._data, 'permanent_interactions',
                                             must_exist=False)
        self._joints = group(self._data, 'joints')
        try:
            self._boundary_conditions = group(self._data, 'boundary_conditions',
                                              must_exist=(self._mode=='w'))
        except Exception as e :
            print('Warning -  group(self._data, boundary_conditions ) : ',  e)
        self._static_data = data(self._data, 'static', 9,
                                 use_compression = self._use_compression)
        self._velocities_data = data(self._data, 'velocities', 8,
                                     use_compression = self._use_compression)
        self._dynamic_data = data(self._data, 'dynamic', 9,
                                  use_compression = self._use_compression)
        self._cf_data = data(self._data, 'cf', 15,
                             use_compression = self._use_compression)
        if self._should_output_domains or 'domain' in self._data:
            self._domain_data = data(self._data, 'domain', 3,
                                     use_compression = self._use_compression)
        self._solv_data = data(self._data, 'solv', 4,
                               use_compression = self._use_compression)
        self._input = group(self._data, 'input')
        self._nslaws_data = group(self._data, 'nslaws')

        if self._shape_filename is None:
            if self._collision_margin:
                self._shape = ShapeCollection(
                    io=self, collision_margin=self._collision_margin)

            else:
                self._shape = ShapeCollection(io=self)
        else:
            if self._collision_margin:
                self._shape = ShapeCollection(
                    io=self._shape_filename,
                    collision_margin=self._collision_margin)
            else:
                self._shape = ShapeCollection(io=self._shape_filename)
        return self

    def __exit__(self, type_, value, traceback):
        self._out.close()

    def apply_gravity(self, body):
        g = constants.g / self._gravity_scale
        weight = [0, 0, - body.scalarMass() * g]
        body.setFExtPtr(weight)

# hdf5 structure

    def shapes(self):
        """
        Shapes : parameterized primitives or user defined
                 (convex set or meshes)
        """
        return self._ref

    def permanent_interactions(self):
        """
        Permanent interactions.
        """
        return self._permanent_interactions

    def static_data(self):
        """
        Coordinates and orientations of static objects.
        """
        return self._static_data

    def dynamic_data(self):
        """
        Coordinates and orientations of dynamic objects.
        """
        return self._dynamic_data

    def velocities_data(self):
        """
        Velocities of dynamic objects
        """
        return self._velocities_data

    def contact_forces_data(self):
        """
        Contact points informations.
        """
        return self._cf_data

    def domains_data(self):
        """
        Contact point domain information.
        """
        return self._domain_data

    def solver_data(self):
        """
        Solver output
        """
        return self._solv_data

    def instances(self):
        """
        Scene objects.
        """
        return self._input

    def nonsmooth_laws(self):
        """
        Non smooth laws between group of contactors.
        """
        return self._nslaws_data

    def joints(self):
        """
        Joints between dynamic objects or between an object and the scenery.
        """
        return self._joints

    def boundary_conditions(self):
        """
        Boundary conditions applied to  dynamic objects
        """
        return self._boundary_conditions

    def importNonSmoothLaw(self, name):
        if self._broadphase is not None:
            nslawClass = getattr(Kernel, self._nslaws_data[name].attrs['type'])
            if nslawClass == Kernel.NewtonImpactFrictionNSL:
                nslaw = nslawClass(float(self._nslaws_data[name].attrs['e']), 0.,
                                   float(self._nslaws_data[name].attrs['mu']), 3)
            elif nslawClass == Kernel.NewtonImpactNSL:
                nslaw = nslawClass(float(self._nslaws_data[name].attrs['e']))
            elif nslawClass == Kernel.RelayNSL:
                nslaw = nslawClass(int(self._nslaws_data[name].attrs['size']),
                                   float(self._nslaws_data[name].attrs['lb']),
                                   float(self._nslaws_data[name].attrs['ub']))
            assert(nslaw)
            self._nslaws[name] = nslaw
            gid1 = int(self._nslaws_data[name].attrs['gid1'])
            gid2 = int(self._nslaws_data[name].attrs['gid2'])
            if gid1 >= 0 and gid2 >= 0:
                if use_proposed:
                    self._broadphase.insertNonSmoothLaw(nslaw, gid1, gid2)
                elif use_original:
                    self._broadphase.insert(nslaw, gid1, gid2)

    def importOccObject(self, name, translation, orientation,
                        velocity, contactors, mass, given_inertia, body_class,
                        shape_class, face_class, edge_class, number=None):

        if mass > 0.:
            if body_class is None:
                body_class = occ.OccBody

            assert (given_inertia is not None)
            inertia = given_inertia

            body = body_class(
                list(translation) + list(orientation), velocity, mass, inertia)

            if number is not None:
                body.setNumber(number)
        else:
            # a static object
            body = None

        ctors_data = set([contactor.data for contactor in contactors])

        ref_shape = {ctor_data: occ.OccContactShape(
            self._shape.get(ctor_data,
                            shape_class, face_class,
                            edge_class, new_instance=True))
                     for ctor_data in ctors_data}

        ref_added = dict()
        for contactor in contactors:

            contact_shape = None
            reference_shape = ref_shape[contactor.data]

            self._keep.append(reference_shape)

            if hasattr(contactor, 'contact_type'):

                if contactor.contact_type == 'Face':
                    contact_shape = \
                                occ.OccContactFace(reference_shape,
                                                   contactor.contact_index)

                elif contactor.contact_type == 'Edge':
                    contact_shape = \
                                occ.OccContactEdge(reference_shape,
                                                   contactor.contact_index)

            if contact_shape is not None:

                if name not in self._occ_contactors:
                    self._occ_contactors[name] = dict()

                self._occ_contactors[name][contactor.instance_name] = \
                                                        contact_shape

                if body is not None:
                    body.addContactShape(contact_shape,
                                         contactor.translation,
                                         contactor.orientation,
                                         contactor.group)
                else:
                    occ.occ_move(contact_shape, list(contactor.translation) +\
                                 list(contactor.orientation))

            if reference_shape not in ref_added:
                if body is not None:
                    body.addShape(reference_shape.shape(),
                                  contactor.translation,
                                  contactor.orientation)
                else:
                    occ.occ_move(reference_shape.shape(),
                                 list(contactor.translation) +\
                                 list(contactor.orientation))
                ref_added[reference_shape] = True

        if body is not None:
            self._set_external_forces(body)

            # add the dynamical system to the non smooth
            # dynamical system
            nsds = self._model.nonSmoothDynamicalSystem()
            nsds.insertDynamicalSystem(body)
            nsds.setName(body, str(name))

    def importBulletObject(self, name, translation, orientation,
                           velocity, contactors, mass, inertia,
                           body_class, shape_class, birth=False,
                           number = None):


        if body_class is None:
            body_class = default_body_class

        if self._broadphase is not None and 'input' in self._data:
            body = None
            if use_proposed and mass == 0:
                # a static object

                cset = SiconosContactorSet()
                csetpos = (translation + orientation)
                for c in contactors:
                    shp = self._shape.get(c.data)
                    pos = list(c.translation) + list(c.orientation)
                    cset.append(SiconosContactor(shp, pos, c.group))
                    print('Adding shape %s to static contactor'%c.data, pos)
                self._broadphase.insertStaticContactorSet(cset, csetpos)

                self._static[name] = {
                    'number': number,
                    'origin': translation,
                    'orientation': orientation,
                    'transform': btTransform(btQuaternion(orientation[1],
                                                          orientation[2],
                                                          orientation[3],
                                                          orientation[0]),
                                             btVector3(translation[0],
                                                       translation[1],
                                                       translation[2])),
                    'shape': shp,
                }

            elif use_original and mass == 0. and use_bullet:
                # a static object
                rbase = btQuaternion(orientation[1],
                                     orientation[2],
                                     orientation[3],
                                     orientation[0])

                tbase = btVector3(translation[0],
                                  translation[1],
                                  translation[2])

                for c in contactors:

                    c_orientation = btQuaternion(c.orientation[1],
                                                 c.orientation[2],
                                                 c.orientation[3],
                                                 c.orientation[0])

                    c_origin = btVector3(c.translation[0],
                                         c.translation[1],
                                         c.translation[2])

                    static_cobj = btCollisionObject()

                    BulletDS.setRelativeTransform(static_cobj,
                                                  tbase,
                                                  rbase,
                                                  c_origin,
                                                  c_orientation)

                    static_cobj.setCollisionFlags(
                        btCollisionObject.CF_STATIC_OBJECT)

                    static_cobj.setCollisionShape(
                        self._shape.get(c.data))

                    self._static[name] = {
                        'number': number,
                        'origin': static_cobj.getWorldTransform().getOrigin(),
                        'orientation': static_cobj.getWorldTransform().getRotation(),
                        'transform': static_cobj.getWorldTransform(),
                        'cobj': static_cobj,
                        }

                    self._broadphase.addStaticObject(static_cobj, int(c.group))

            elif use_proposed:
                # a proposed-API moving object

                if inertia is not None:
                    if np.shape(inertia) == (3,):
                        inertia = np.diag(inertia)
                    elif np.shape(inertia) != (3,3):
                        print('Wrong shape of inertia')
                    have_inertia = True
                else:
                    have_inertia = False

                body = body_class(translation + orientation,
                                  velocity,
                                  mass, inertia)
                if have_inertia:
                    body.setUseContactorInertia(False)

                self_collide = self._input[name].get('allow_self_collide',None)
                if self_collide is not None:
                    body.setAllowSelfCollide(not not self_collide)

                cset = SiconosContactorSet()
                for c in contactors:
                    shp = self._shape.get(c.data)
                    pos = list(c.translation) + list(c.orientation)
                    cset.append(SiconosContactor(shp, pos, c.group))

                body.setContactors(cset)

            elif use_original and use_bullet:
                # a Bullet moving object
                bws = BulletWeightedShape(
                    self._shape.get(contactors[0].data), mass)

                if inertia is not None:

                    if np.shape(inertia) == (3,):
                        bws.setInertia(inertia[0], inertia[1], inertia[2])
                    elif (np.shape(inertia) == (3, 3)):
                        bws.setInertia(inertia)
                    else:
                        print('Wrong shape of inertia')

                body = body_class(bws,
                                  translation + orientation,
                                  velocity,
                                  contactors[0].translation,
                                  contactors[0].orientation,
                                  contactors[0].group)

                #body.setNullifyFGyr(True)
                for contactor in contactors[1:]:
                    shape_id = self._shapeid[contactor.data]

                    body.addCollisionShape(self._shape.get(contactor.data),
                                           contactor.translation,
                                           contactor.orientation,
                                           contactor.group)

            if body:
                # set id number
                if number is not None:
                    body.setNumber(number)

                # set external forces
                self._set_external_forces(body)

                # add the dynamical system to the non smooth
                # dynamical system
                if birth:
                    nsds = self._model.nonSmoothDynamicalSystem()
                    if use_proposed:
                        nsds.insertDynamicalSystem(body)
                        self._model.simulation().prepareIntegratorForDS(
                            self._osi, body, self._model,
                            self._model.simulation().nextTime())
                        self._model.simulation().initialize(self._model, False)
                    elif use_original:
                        self._broadphase.addDynamicObject(
                            body,
                            self._model.simulation(),
                            self._osi)
                    nsds.setName(body, str(name))
                else:
                    nsds = self._model.nonSmoothDynamicalSystem()
                    nsds.insertDynamicalSystem(body)
                    nsds.setName(body, str(name))

    def importJoint(self, name):
        if self._broadphase is not None:
            nsds = self._model.nonSmoothDynamicalSystem()
            topo = nsds.topology()

            joint_type = self.joints()[name].attrs['type']
            joint_class = getattr(joints, joint_type)
            absolute = self.joints()[name].attrs.get('absolute', None)
            absolute = [[True if absolute else False], []][absolute is None]
            allow_self_collide = self.joints()[name].attrs.get(
                'allow_self_collide',None)
            stops = self.joints()[name].attrs.get('stops',None)
            nslaws = self.joints()[name].attrs.get('nslaws',None)
            friction = self.joints()[name].attrs.get('friction',None)

            ds1_name = self.joints()[name].attrs['object1']
            ds1 = topo.getDynamicalSystem(ds1_name)
            ds2 = None

            if 'object2' in self.joints()[name].attrs:
                ds2_name = self.joints()[name].attrs['object2']
                ds2 = topo.getDynamicalSystem(ds2_name)
                try:
                    joint = joint_class(ds1, ds2,
                                        self.joints()[name].attrs['pivot_point'],
                                        self.joints()[name].attrs['axis'],
                                        *absolute)
                except NotImplementedError:
                    try:
                        joint = joint_class(ds1, ds2,
                                            self.joints()[name].attrs['pivot_point'],
                                            *absolute)
                    except NotImplementedError:
                        joint = joint_class(ds1, ds2, *absolute)

            else:
                try:
                    joint = joint_class(ds1,
                                        self.joints()[name].attrs['pivot_point'],
                                        self.joints()[name].attrs['axis'],
                                        *absolute)
                except NotImplementedError:
                    try:
                        joint = joint_class(ds1,
                                            self.joints()[name].attrs['pivot_point'],
                                            *absolute)
                    except NotImplementedError:
                        joint = joint_class(ds1, *absolute)

            if allow_self_collide is not None:
                joint.setAllowSelfCollide(not not allow_self_collide)
            joint_nslaw = EqualityConditionNSL(joint.numberOfConstraints())
            joint_inter = Interaction(joint_nslaw, joint)
            self._model.nonSmoothDynamicalSystem().\
                link(joint_inter, ds1, ds2)
            nsds.setName(joint_inter, str(name))

            # Add a e=0 joint by default, otherwise user can specify
            # the impact law by name or a list of names for each axis.
            if stops is not None:
                assert np.shape(stops)[1] == 3, 'Joint stops shape must be (?,3)'
                if nslaws is None:
                    nslaws = [Kernel.NewtonImpactNSL(0.0)]*np.shape(stops)[0]
                elif isinstance(nslaws,str):
                    nslaws = [self._nslaws[nslaws]]*np.shape(stops)[0]
                else:
                    assert(np.shape(nslaws)[0]==np.shape(stops)[0])
                    nslaws = [self._nslaws[nsl] for nsl in nslaws]
                for n, (nsl, (axis, pos, dir)) in enumerate(zip(nslaws,stops)):
                    # "bool()" is needed because type of dir is
                    # numpy.bool_, which SWIG doesn't handle well.
                    stop = joints.JointStopR(joint, pos, bool(dir<0), int(axis))
                    stop_inter = Interaction(nsl, stop)
                    self._model.nonSmoothDynamicalSystem().\
                        link(stop_inter, ds1, ds2)
                    nsds.setName(stop_inter, '%s_stop%d'%(str(name),n))

            # The per-axis friction NSL, can be ''
            if friction is not None:
                if isinstance(friction,basestring):
                    friction = [friction]
                else:
                    assert hasattr(friction, '__iter__')
                for ax,fr_nslaw in enumerate(friction):
                    if fr_nslaw == '':  # no None in hdf5, use empty string
                        continue        # instead for no NSL on an axis
                    nslaw = self._nslaws[fr_nslaw]
                    fr = joints.JointFrictionR(joint, [ax])
                    fr_inter = Interaction(nslaw, fr)
                    self._model.nonSmoothDynamicalSystem().\
                        link(fr_inter, ds1, ds2)
                    nsds.setName(fr_inter, '%s_friction%d'%(str(name),ax))

    def importBoundaryConditions(self, name):
        if self._broadphase is not None:
            topo = self._model.nonSmoothDynamicalSystem().\
                topology()

            bc_type = self.boundary_conditions()[name].attrs['type']
            bc_class = getattr(Kernel,bc_type)

            print('name = ', name)
            print('object1')

            ds1_name = self.boundary_conditions()[name].attrs['object1']
            ds1 = topo.getDynamicalSystem(ds1_name)


            if ( bc_type == 'HarmonicBC') :
                bc = bc_class(self.boundary_conditions()[name].attrs['indices'],
                              self.boundary_conditions()[name].attrs['a'],
                              self.boundary_conditions()[name].attrs['b'],
                              self.boundary_conditions()[name].attrs['omega'],
                              self.boundary_conditions()[name].attrs['phi'])

            elif ( bc_type == 'FixedBC' ):
                bc = bc_class(self.boundary_conditions()[name].attrs['indices'])

            elif ( bc_type == 'BoundaryCondition' ):
                bc = bc_class(self.boundary_conditions()[name].attrs['indices'],
                              self.boundary_conditions()[name].attrs['v'])

            # set bc to the ds1

            ds1.setBoundaryConditions(bc);

            #joint_inter = Interaction(joint_nslaw, joint)
            #    self._model.nonSmoothDynamicalSystem().\
            #        link(joint_inter, ds1)

    def importPermanentInteraction(self, name):
        """
        """
        if (self._broadphase is not None and 'input' in self._data
              and self.permanent_interactions() is not None):
            topo = self._model.nonSmoothDynamicalSystem().\
                topology()

            pinter = self.permanent_interactions()[name]
            body1_name = pinter.attrs['body1_name']
            body2_name = pinter.attrs['body2_name']

            ds1 = occ.cast_OccBody(topo.getDynamicalSystem(body1_name))
            try:
                ds2 = occ.cast_OccBody(topo.getDynamicalSystem(body2_name))
            except:
                ds2 = None

            contactor1_name = pinter.attrs['contactor1_name']
            contactor2_name = pinter.attrs['contactor2_name']

            distance_calculator = pinter.attrs['distance_calculator']
            offset = pinter.attrs['offset']

            body1 = self._input[body1_name]
            body2 = self._input[body2_name]

            ctr1 = body1[contactor1_name]
            ctr2 = body2[contactor2_name]

            cg1 = int(ctr1.attrs['group'])
            cg2 = int(ctr2.attrs['group'])
            nslaw = self._broadphase.nslaw(cg1, cg2)

            cocs1 = self._occ_contactors[body1_name][contactor1_name]
            cocs2 = self._occ_contactors[body2_name][contactor2_name]

            # else:
            #     topods2 = self._shape.get(ctr2.attrs['name'])
            #     self._keep.append(topods2)
            #     ocs2 = occ.OccContactShape(topods2)

            #     index2 = int(ctr2.attrs['contact_index'])

            #     ctact_t2 = ctr2.attrs['type']

            #     ctactbuild = {'Face': occ.OccContactFace,
            #                   'Edge': occ.OccContactEdge}

            #     cocs2 = ctactbuild[ctact_t2](ocs2, index2)

            cp1 = occ.ContactPoint(cocs1)
            cp2 = occ.ContactPoint(cocs2)

            real_dist_calc = {'cadmbtb': occ.CadmbtbDistanceType,
                              'occ': occ.OccDistanceType}

            relation = occ.OccR(cp1, cp2,
                                real_dist_calc[distance_calculator]())

            relation.setOffset(offset)

            inter = Interaction(nslaw, relation)

            if ds2 is not None:
                self._model.nonSmoothDynamicalSystem().link(inter, ds1, ds2)
            else:
                self._model.nonSmoothDynamicalSystem().link(inter, ds1)

            # keep pointers
            self._keep.append([cocs1, cocs2, cp1,
                               cp2, relation])

    def importScene(self, time, body_class, shape_class, face_class,
                    edge_class):
        """
        From the specification given in the hdf5 file with the help of
        add* functions, import into the broadphase object:
          - the static objects
          - the dynamic objects
          - the joints
          - the nonsmooth laws
        that have a specified time of birth <= current time.
        """

        # Ensure we count up from zero for implicit DS numbering
        DynamicalSystem.resetCount(0)

        for shape_name in self._ref:
            self._shapeid[shape_name] = self._ref[shape_name].attrs['id']
            self._number_of_shapes += 1

        # import dynamical systems
        if self._broadphase is not None and 'input' in self._data:

            dpos_data = self.dynamic_data()
            if dpos_data is not None and len(dpos_data) > 0:

                max_time = max(dpos_data[:, 0])
                id_last = np.where(
                    abs(dpos_data[:, 0] - max_time) < 1e-9)[0]

            else:
                # should not be used
                max_time = None
                id_last = None

            for (name, obj) in sorted(self._input.items(),
                                      key=lambda x: x[0]):

                mass = obj.attrs['mass']
                time_of_birth = obj.attrs['time_of_birth']

                if time_of_birth >= time:
                    #
                    # in the future
                    #
                    bisect.insort_left(self._scheduled_births, time_of_birth)
                    if time_of_birth in self._births:
                        self._births[time_of_birth].append((name, obj))
                    else:
                        self._births[time_of_birth] = [(name, obj)]
                else:
                    #
                    # this is for now
                    #
                    # cold restart if output previously done
                    if mass > 0 and dpos_data is not None and len(dpos_data) > 0:
                        print ('Import  dynamic object name ', name, 'from current state')
                        print ('  number of imported object ', obj.attrs['id'])

                        id_last_inst = np.where(
                            dpos_data[id_last, 1] ==
                            self.instances()[name].attrs['id'])[0]
                        xpos = dpos_data[id_last[id_last_inst[0]], :]
                        translation = (xpos[2], xpos[3], xpos[4])
                        orientation = (xpos[5], xpos[6], xpos[7], xpos[8])

                        velocities = self.velocities_data()
                        id_vlast = np.where(
                            abs(velocities[:, 0] - max_time) < 1e-9)[0]

                        id_vlast_inst = np.where(
                            velocities[id_vlast, 1] ==
                            self.instances()[name].attrs['id'])[0]
                        xvel = velocities[id_vlast[id_vlast_inst[0]], :]
                        velocity = (xvel[2], xvel[3], xvel[4], xvel[5], xvel[6], xvel[7])

                    # start from initial conditions
                    else:
                        print ('Import  dynamic or static object number ', obj.attrs['id'], 'from initial state')
                        print ('                object name   ', name)
                        translation = obj.attrs['translation']
                        orientation = obj.attrs['orientation']
                        velocity = obj.attrs['velocity']

                    # bodyframe center of mass
                    center_of_mass = None

                    # bodyframe center of mass
                    # check for compatibility
                    if 'center_of_mass' in obj.attrs:
                        center_of_mass = \
                                obj.attrs['center_of_mass'].astype(float)
                    else:
                        center_of_mass = [0, 0, 0]

                    input_ctrs = [ctr for _n_, ctr in obj.items()]

                    contactors = []
                    occ_type = False
                    for ctr in input_ctrs:
                        if 'type' in ctr.attrs:
                            # occ contact
                            occ_type = True
                            contactors.append(
                                Contactor(
                                    instance_name=ctr.attrs['instance_name'],
                                    shape_data=ctr.attrs['name'],
                                    collision_group=ctr.attrs['group'].astype(int),
                                    contact_type=ctr.attrs['type'],
                                    contact_index=ctr.attrs['contact_index'].astype(int),
                                    relative_translation=np.subtract(ctr.attrs['translation'].astype(float), center_of_mass),
                                    relative_orientation=ctr.attrs['orientation'].astype(float)))
                        elif 'group' in ctr.attrs:
                            # bullet contact
                            assert not occ_type
                            contactors.append(
                                Contactor(
                                    instance_name=ctr.attrs['instance_name'],
                                    shape_data=ctr.attrs['name'],
                                    collision_group=ctr.attrs['group'].astype(int),
                                    relative_translation=np.subtract(ctr.attrs['translation'].astype(float), center_of_mass),
                                    relative_orientation=ctr.attrs['orientation'].astype(float)))
                        else:
                            # occ shape
                            occ_type = True
                            # fix: not only contactors here
                            contactors.append(
                                Volume(
                                    instance_name=ctr.attrs['instance_name'],
                                    shape_data=ctr.attrs['name'],
                                    parameters=pickle.loads(ctr.attrs['parameters']),
                                    relative_translation=np.subtract(ctr.attrs['translation'].astype(float), center_of_mass),
                                    relative_orientation=ctr.attrs['orientation'].astype(float)))

                    if 'inertia' in obj.attrs:
                        inertia = obj.attrs['inertia']
                    else:
                        inertia = None

                    if occ_type:
                        # Occ object
                        self.importOccObject(
                            name, floatv(translation), floatv(orientation),
                            floatv(velocity), contactors, float(mass),
                            inertia, body_class, shape_class, face_class,
                            edge_class,
                            number = self.instances()[name].attrs['id'])
                    else:
                        # Bullet object
                        self.importBulletObject(
                            name, floatv(translation), floatv(orientation),
                            floatv(velocity), contactors, float(mass),
                            inertia, body_class, shape_class,
                            number = self.instances()[name].attrs['id'])
            # import nslaws
            # note: no time of birth for nslaws and joints
            for name in self._nslaws_data:
                self.importNonSmoothLaw(name)

            for name in self.joints():
                self.importJoint(name)

            for name in self.boundary_conditions():
                self.importBoundaryConditions(name)

            for name in self.permanent_interactions():
                self.importPermanentInteraction(name)

    def currentTime(self):
        if self._initializing:
            return self._model.simulation().startingTime()
        else:
            return self._model.simulation().nextTime()

    def importBirths(self, body_class=None, shape_class=None,
                     face_class=None, edge_class=None,):
        """
        Import new objects in the broadphase.
        """
        time = self.currentTime()

        ind_time = bisect.bisect_left(self._scheduled_births, time)

        current_times_of_births = set(self._scheduled_births[:ind_time])
        self._scheduled_births = self._scheduled_births[ind_time:]

        #print (time, current_times_of_births)
        for time_of_birth in current_times_of_births:
            #print( "time_of_birth", time_of_birth)
            for (name, obj) in self._births[time_of_birth]:
                #print(name,obj)
                translation = obj.attrs['translation']
                orientation = obj.attrs['orientation']
                velocity = obj.attrs['velocity']

                input_ctrs = [ctr for _n_, ctr in obj.items()]
                mass = obj.attrs['mass']

                contactors = [Contactor(
                    instance_name=ctr.attrs['instance_name'],
                    shape_data=ctr.attrs['name'],
                    collision_group=int(ctr.attrs['group']),
                    relative_translation=floatv(ctr.attrs['translation']),
                    relative_orientation=floatv(ctr.attrs['orientation']))
                              for ctr in input_ctrs]

                if 'inertia' in obj.attrs:
                    inertia = obj.attrs['inertia']
                else:
                    inertia = None

                if True in ('type' in self.shapes()[ctr.attrs['name']].attrs
                            and self.shapes()[ctr.attrs['name']].attrs['type']
                            in ['brep', 'step']
                            for ctr in input_ctrs):
                    # Occ object
                    self.importOccObject(
                        name, floatv(translation), floatv(orientation),
                        floatv(
                            velocity), contactors, float(mass),
                        inertia, body_class, shape_class, face_class,
                        edge_class, number = self.instances()[name].attrs['id'])
                else:
                    # Bullet object
                    self.importBulletObject(
                        name, floatv(translation), floatv(orientation),
                        floatv(velocity), contactors, float(mass),
                        inertia, body_class, shape_class, birth=True,
                        number = self.instances()[name].attrs['id'])

    def outputStaticObjects(self):
        """
        Outputs translations and orientations of static objects
        """
        time = self.currentTime()
        p = 0
        self._static_data.resize(len(self._static), 0)

        for static in self._static.values():
            print('output static object', static['number'])
            translation = static['transform'].getOrigin()
            rotation = static['transform'].getRotation()
            self._static_data[p, :] = \
                [time,
                 static['number'],
                 translation.x(),
                 translation.y(),
                 translation.z(),
                 rotation.w(),
                 rotation.x(),
                 rotation.y(),
                 rotation.z()]
            p += 1

    def outputDynamicObjects(self, initial=False):
        """
        Outputs translations and orientations of dynamic objects.
        """

        current_line = self._dynamic_data.shape[0]

        time = self.currentTime()

        positions = self._io.positions(self._model)

        if positions is not None:

            self._dynamic_data.resize(current_line + positions.shape[0], 0)

            times = np.empty((positions.shape[0], 1))
            times.fill(time)

            self._dynamic_data[current_line:, :] = np.concatenate((times,
                                                                   positions),
                                                                  axis=1)

    def outputVelocities(self):
        """
        Output velocities of dynamic objects
        """

        current_line = self._dynamic_data.shape[0]

        time = self.currentTime()

        velocities = self._io.velocities(self._model)

        if velocities is not None:

            self._velocities_data.resize(current_line + velocities.shape[0], 0)

            times = np.empty((velocities.shape[0], 1))
            times.fill(time)

            self._velocities_data[current_line:, :] = np.concatenate((times,
                                                                      velocities),
                                                                     axis=1)

    def outputContactForces(self):
        """
        Outputs contact forces
        _contact_index_set default value is 1.
        """
        if self._model.nonSmoothDynamicalSystem().\
                topology().indexSetsSize() > 1:
            time = self.currentTime()
            contact_points = self._io.contactPoints(self._model,
                                                    self._contact_index_set)

            if contact_points is not None:

                current_line = self._cf_data.shape[0]
                self._cf_data.resize(current_line + contact_points.shape[0], 0)
                times = np.empty((contact_points.shape[0], 1))
                times.fill(time)

                self._cf_data[current_line:, :] = \
                    np.concatenate((times,
                                    contact_points),
                                   axis=1)

    def outputDomains(self):
        """
        Outputs domains of contact points
        """
        if self._model.nonSmoothDynamicalSystem().\
                topology().indexSetsSize() > 1:
            time = self.currentTime()
            domains = self._io.domains(self._model)

            if domains is not None:

                current_line = self._domain_data.shape[0]
                self._domain_data.resize(current_line + domains.shape[0], 0)
                times = np.empty((domains.shape[0], 1))
                times.fill(time)

                self._domain_data[current_line:, :] = \
                    np.concatenate((times, domains), axis=1)

    def outputSolverInfos(self):
        """
        Outputs solver #iterations & precision reached
        """

        time = self.currentTime()
        so = self._model.simulation().oneStepNSProblem(0).\
            numericsSolverOptions()

        current_line = self._solv_data.shape[0]
        self._solv_data.resize(current_line + 1, 0)
        if so.solverId == Numerics.SICONOS_GENERIC_MECHANICAL_NSGS:
            iterations = so.iparam[3]
            precision = so.dparam[2]
            local_precision = so.dparam[3]
        elif so.solverId == Numerics.SICONOS_FRICTION_3D_NSGS:
            iterations = so.iparam[Numerics.SICONOS_IPARAM_ITER_DONE]
            precision = so.dparam[Numerics.SICONOS_DPARAM_RESIDU]
            local_precision = 0.
        # maybe wrong for others
        else:
            iterations = so.iparam[Numerics.SICONOS_IPARAM_ITER_DONE]
            precision = so.dparam[Numerics.SICONOS_DPARAM_RESIDU]
            local_precision = so.dparam[2]

        self._solv_data[current_line, :] = [time, iterations, precision,
                                            local_precision]

    def printSolverInfos(self):
        """
        Outputs solver #iterations & precision reached
        """
        time = self.currentTime()
        so = self._model.simulation().oneStepNSProblem(0).\
            numericsSolverOptions()
        if so.solverId == Numerics.SICONOS_GENERIC_MECHANICAL_NSGS:
            iterations = so.iparam[3]
            precision = so.dparam[2]
            local_precision = so.dparam[3]
        elif so.solverId == Numerics.SICONOS_FRICTION_3D_NSGS:
            iterations = so.iparam[Numerics.SICONOS_IPARAM_ITER_DONE]
            precision = so.dparam[Numerics.SICONOS_DPARAM_RESIDU]
            local_precision = 0.
        # maybe wrong for others
        else:
            iterations = so.iparam[Numerics.SICONOS_IPARAM_ITER_DONE]
            precision = so.dparam[Numerics.SICONOS_DPARAM_RESIDU]
            local_precision = so.dparam[2]


        print('SolverInfos at time :', time,
              'iterations= ', iterations,
              'precision=', precision,
              'local_precision=', )

    def addMeshFromString(self, name, shape_data, scale=None,
                          insideMargin=None, outsideMargin=None):
        """
        Add a mesh shape from a string.
        Accepted format : mesh encoded in VTK .vtp format
        """

        import vtk

        if name not in self._ref:

            shape = self._ref.create_dataset(name, (1,),
                                             dtype=h5py.new_vlen(str))
            shape[:] = shape_data
            shape.attrs['id'] = self._number_of_shapes
            shape.attrs['type'] = 'vtp'
            if scale is not None:
                shape.attrs['scale'] = scale
            if insideMargin is not None:
                shape.attrs['insideMargin'] = insideMargin
            if outsideMargin is not None:
                shape.attrs['outsideMargin'] = outsideMargin
            self._shapeid[name] = shape.attrs['id']
            self._number_of_shapes += 1

    def addMeshFromFile(self, name, filename, scale=None,
                        insideMargin=None, outsideMargin=None):
        """
        Add a mesh shape from a file.
        Accepted format : .stl or mesh encoded in VTK .vtp format
        """

        import vtk

        if filename[0] != os.path.sep:
            filename = os.path.join(os.path.split(os.path.abspath(sys.argv[0]))[0],
                                    filename)
        if name not in self._ref:

            if os.path.splitext(filename)[-1][1:] == 'stl':
                reader = vtk.vtkSTLReader()
                reader.SetFileName(filename)
                reader.Update()

                if reader.GetErrorCode() != 0:
                    print('vtkSTLReader error', reader.GetErrorCode())
                    sys.exit(1)

                with tmpfile() as tmpf:
                    writer = vtk.vtkXMLPolyDataWriter()
                    writer.SetInputData(reader.GetOutput())
                    writer.SetFileName(tmpf[1])
                    writer.Write()

                    shape_data = str_of_file(tmpf[1])

            else:
                assert os.path.splitext(filename)[-1][1:] == 'vtp'
                shape_data = str_of_file(filename)

            self.addMeshFromString(name, shape_data, scale=scale,
                                   insideMargin=insideMargin,
                                   outsideMargin=outsideMargin)

    def addHeightMap(self, name, heightmap, rectangle,
                     insideMargin=None, outsideMargin=None):
        """
        Add a heightmap represented as a SiconosMatrix
        """
        assert(heightmap.shape[0] >= 2)
        assert(heightmap.shape[1] >= 2)
        if name not in self._ref:
            shape = self._ref.create_dataset(name, data=heightmap)
            shape.attrs['id'] = self._number_of_shapes
            shape.attrs['type'] = 'heightmap'

            # measurements of the heightfield, i.e. length of sides of
            # the rectangle where heightmap will be placed -- height
            # is represented by heightmap values
            assert(len(rectangle)==2)
            shape.attrs['rect'] = rectangle # tuple (length x, length y)

            if insideMargin is not None:
                shape.attrs['insideMargin'] = insideMargin
            if outsideMargin is not None:
                shape.attrs['outsideMargin'] = outsideMargin
            self._shapeid[name] = shape.attrs['id']
            self._number_of_shapes += 1

    def addBRepFromString(self, name, shape_data):
        """
        Add a brep contained in a string.
        """
        if name not in self._ref:
            shape = self._ref.create_dataset(name, (1,),
                                             dtype=h5py.new_vlen(str))
            if type(shape_data) == str:
                # raw str
                shape[:] = shape_data
            else:
                # __getstate__ as with pythonocc
                shape[:] = shape_data[0]
                shape.attrs['occ_indx'] = shape_data[1]

            shape.attrs['id'] = self._number_of_shapes
            shape.attrs['type'] = 'brep'

            self._shapeid[name] = shape.attrs['id']
            self._number_of_shapes += 1

    def addOccShape(self, name, occ_shape):
        """
        Add an OpenCascade TopoDS_Shape.
        """

        if name not in self._ref:

            from OCC.STEPControl import STEPControl_Writer, STEPControl_AsIs

            # step format is used for the storage.
            step_writer = STEPControl_Writer()

            step_writer.Transfer(occ_shape, STEPControl_AsIs)

            shape_data = None

            with tmpfile() as tmpf:

                status = step_writer.Write(tmpf[1])

                tmpf[0].flush()
                shape_data = str_of_file(tmpf[1])

                shape = self._ref.create_dataset(name, (1,),
                                                 dtype=h5py.new_vlen(str))
                shape[:] = shape_data
                shape.attrs['id'] = self._number_of_shapes
                shape.attrs['type'] = 'step'
                self._shapeid[name] = shape.attrs['id']
                self._number_of_shapes += 1

    def addShapeDataFromFile(self, name, filename):
        """
        Add shape data from a file.
        """
        if name not in self._ref:
            shape = self._ref.create_dataset(name, (1,),
                                             dtype=h5py.new_vlen(str))
            shape[:] = str_of_file(filename)
            shape.attrs['id'] = self._number_of_shapes
            try:
                shape.attrs['type'] = os.path.splitext(filename)[1][1:]
            except:
                shape.attrs['type'] = 'unknown'

            self._shapeid[name] = shape.attrs['id']
            self._number_of_shapes += 1

    def addInteraction(self, name, body1_name, contactor1_name,
                       body2_name, contactor2_name,
                       distance_calculator='cadmbtb',
                       offset=0.0001):
        """
        Add permanent interactions between two objects contactors.
        """
        if name not in self.permanent_interactions():
            pinter = self.permanent_interactions().\
                      create_dataset(name, (1,),
                                     dtype=h5py.new_vlen(str))
            pinter.attrs['id'] = self._number_of_permanent_interactions
            pinter.attrs['type'] = 'permanent_interaction'
            pinter.attrs['body1_name'] = body1_name
            pinter.attrs['body2_name'] = body2_name
            pinter.attrs['contactor1_name'] = contactor1_name
            pinter.attrs['contactor2_name'] = contactor2_name
            pinter.attrs['distance_calculator'] = distance_calculator
            pinter.attrs['offset'] = offset

            self._pinterid[name] = pinter.attrs['id']
            self._number_of_permanent_interactions += 1

    def addConvexShape(self, name, points,
                       insideMargin=None, outsideMargin=None):
        """
        Add a convex shape defined by a list of points.
        """
        if name not in self._ref:
            shape=self._ref.create_dataset(name,
                                             (np.shape(points)[0],
                                              np.shape(points)[1]))
            if insideMargin is not None:
                shape.attrs['insideMargin'] = insideMargin
            if outsideMargin is not None:
                shape.attrs['outsideMargin'] = outsideMargin
            shape[:]=points[:]
            shape.attrs['type']='convex'
            shape.attrs['id']=self._number_of_shapes
            self._shapeid[name]=shape.attrs['id']
            self._number_of_shapes += 1

    def addPrimitiveShape(self, name, primitive, params,
                          insideMargin=None, outsideMargin=None):
        """
        Add a primitive shape.
        """
        if name not in self._ref:
            shape=self._ref.create_dataset(name, (1, len(params)))
            shape.attrs['id']=self._number_of_shapes
            shape.attrs['type']='primitive'
            shape.attrs['primitive']=primitive
            if insideMargin is not None:
                shape.attrs['insideMargin'] = insideMargin
            if outsideMargin is not None:
                shape.attrs['outsideMargin'] = outsideMargin
            shape[:]=params
            self._shapeid[name]=shape.attrs['id']
            self._number_of_shapes += 1

    def addObject(self, name, shapes,
                  translation,
                  orientation=[1, 0, 0, 0],
                  velocity=[0, 0, 0, 0, 0, 0],
                  mass=0, center_of_mass=[0, 0, 0],
                  inertia=None, time_of_birth=-1,
                  allow_self_collide=False):
        """Add an object with associated shapes as a list of Volume or
        Contactor objects. Contact detection and processing is
        defined by the Contactor objects. The Volume objects are used for
        the computation of inertia and center of mass if not provided.

        Each Contactor and Volume object may have a relative
        translation and a relative orientation expressed in the bodyframe
        coordinates.

        Parameters
        ----------
        name: string
            The name of the object.

        shapes: iterable
            The list of associated Contactor or Volume objects.

        translation: array_like of length 3
            Initial translation of the object (mandatory)

        velocity: array_like of length 6
            Initial velocity of the object.
            The components are those of the translation velocity along
            x, y and z axis and the rotation velocity around x, y and
            z axis.  The default velocity is [0, 0, 0, 0, 0, 0].

        mass: float
            The mass of the object, if it is zero the object is defined as
            a static object involved only in contact detection.
            The default value is zero.

        center_of_mass: array_like of length 3
            The position of the center of mass expressed in the body frame
            coordinates.

        inertia: array_like of length 3 or 3x3 matrix.
            The principal moments of inertia (array of length 3) or
            a full 3x3 inertia matrix

        """
        # print(arguments())
        if len(orientation) == 2:
            # axis + angle
            axis=orientation[0]
            assert len(axis) == 3
            angle=orientation[1]
            assert type(angle) is float
            n=sin(angle / 2.) / np.linalg.norm(axis)

            ori=[cos(angle / 2.), axis[0] * n, axis[1] * n, axis[2] * n]
        else:
            assert(len(orientation)==4)
            # a given quaternion
            ori=orientation

        assert(len(translation)==3)

        if name not in self._input:

            com_translation = [0., 0., 0.]

            if inertia is None and mass > 0.:
                volumes = filter(lambda s: isinstance(s, Volume),
                                 shapes)
                if len(volumes) > 0:
                    # a computed inertia and center of mass
                    # occ only
                    from OCC.GProp import GProp_GProps
                    from OCC.BRepGProp import brepgprop_VolumeProperties
                    from OCC.gp import gp_Ax1, gp_Dir
                    volumes = filter(lambda s: isinstance(s, Volume),
                                     shapes)

                    props = GProp_GProps()

                    for volume in volumes:

                        iprops = GProp_GProps()
                        iishape = self._shape.get(volume.data,
                                                  shape_class=None, face_class=None,
                                                  edge_class=None, new_instance=True)
                        ishape = occ.OccContactShape(iishape).data()

                        # the shape relative displacement
                        occ.occ_move(ishape, list(volume.translation) +\
                                     list(volume.orientation))

                        brepgprop_VolumeProperties(iishape, iprops)

                        if volume.parameters is not None and \
                           hasattr(volume.parameters, 'density'):
                            density = volume.parameters.density
                        else:
                            density = 1

                        props.Add(iprops, density)

                    assert (props.Mass() > 0.)
                    global_density = mass / props.Mass()
                    computed_com = props.CentreOfMass()
                    I1 = global_density * props.MomentOfInertia(
                        gp_Ax1(computed_com, gp_Dir(1, 0, 0)))
                    I2 = global_density * props.MomentOfInertia(
                        gp_Ax1(computed_com, gp_Dir(0, 1, 0)))
                    I3 = global_density * props.MomentOfInertia(
                        gp_Ax1(computed_com, gp_Dir(0, 0, 1)))

                    inertia = [I1, I2, I3]
                    center_of_mass = np.array([computed_com.Coord(1),
                                               computed_com.Coord(2),
                                               computed_com.Coord(3)])

                    print('computed inertia:', I1, I2, I3)
                    print('computed center of mass:',
                          computed_com.Coord(1),
                          computed_com.Coord(2),
                          computed_com.Coord(3))


            obj=group(self._input, name)

            obj.attrs['time_of_birth']=time_of_birth

            obj.attrs['mass']=mass
            obj.attrs['translation']=translation
            obj.attrs['orientation']=ori
            obj.attrs['velocity']=velocity
            obj.attrs['center_of_mass']=center_of_mass

            if inertia is not None:
                obj.attrs['inertia']=inertia
            if allow_self_collide is not None:
                obj.attrs['allow_self_collide']=allow_self_collide

            contactors = shapes

            for num, ctor in enumerate(contactors):

                if ctor.instance_name is not None:
                    # a specified name
                    instance_name = ctor.instance_name
                else:
                    # the default name for contactor
                    instance_name = '{0}-{1}'.format(ctor.data, num)

                dat = data(obj, instance_name, 0,
                           use_compression=self._use_compression)

                dat.attrs['instance_name'] = instance_name
                dat.attrs['name'] = ctor.data
                if hasattr(ctor, 'group'):
                    dat.attrs['group'] = ctor.group

                if hasattr(ctor, 'parameters') and \
                        ctor.parameters is not None:
                    dat.attrs['parameters'] = pickle.dumps(ctor.parameters)

                if hasattr(ctor, 'contact_type') and \
                   ctor.contact_type is not None:
                    dat.attrs['type'] = ctor.contact_type

                if hasattr(ctor, 'contact_index') and \
                   ctor.contact_index is not None:
                    dat.attrs['contact_index'] = ctor.contact_index

                dat.attrs['translation'] = ctor.translation
                dat.attrs['orientation'] = ctor.orientation

            if mass == 0:
                obj.attrs['id']=- (self._number_of_static_objects + 1)
                self._number_of_static_objects += 1

            else:
                obj.attrs['id']=(self._number_of_dynamic_objects + 1)
                self._number_of_dynamic_objects += 1

            return obj.attrs['id']

    def addNewtonImpactFrictionNSL(self, name, mu, e=0, collision_group1=0,
                                   collision_group2=0):
        """
        Add a nonsmooth law for contact between 2 groups.
        Only NewtonImpactFrictionNSL are supported.
        name is an user identifiant and must be unique,
        mu is the coefficient of friction,
        e is the coefficient of restitution on the contact normal,
        gid1 and gid2 define the group identifiants.

        """
        if name not in self._nslaws_data:
            nslaw=self._nslaws_data.create_dataset(name, (0,))
            nslaw.attrs['type']='NewtonImpactFrictionNSL'
            nslaw.attrs['mu']=mu
            nslaw.attrs['e']=e
            nslaw.attrs['gid1']=collision_group1
            nslaw.attrs['gid2']=collision_group2

    # Note, default groups are -1 here, indicating not to add them to
    # the nslaw lookup table for contacts, since 1D impacts are
    # useless in this case.  They are however useful for joint stops.
    def addNewtonImpactNSL(self, name, e=0, collision_group1=-1,
                           collision_group2=-1):
        """
        Add a nonsmooth law for contact between 2 groups.
        Only NewtonImpactNSL are supported.
        name is a user identifier and must be unique,
        e is the coefficient of restitution on the contact normal,
        gid1 and gid2 define the group identifiers.

        As opposed to addNewtonImpactFrictionNSL, the default groups are
        -1, making the NSL unassociated with point contacts.  It can
        by used for joint stops however.
        """
        if name not in self._nslaws_data:
            nslaw=self._nslaws_data.create_dataset(name, (0,))
            nslaw.attrs['type']='NewtonImpactNSL'
            nslaw.attrs['e']=e
            nslaw.attrs['gid1']=collision_group1
            nslaw.attrs['gid2']=collision_group2

    # Note, default groups are -1 here, indicating not to add them to
    # the nslaw lookup table for contacts, since 1D impacts are
    # useless in this case.  They are however useful for joint friction.
    def addRelayNSL(self, name, lb, ub, size=1, collision_group1=-1,
                    collision_group2=-1):
        """
        Add a nonsmooth law for contact between 2 groups.
        Only NewtonImpactNSL are supported.
        name is a user identifier and must be unique,
        e is the coefficient of restitution on the contact normal,
        gid1 and gid2 define the group identifiers.

        As opposed to addNewtonImpactFrictionNSL, the default groups are
        -1, making the NSL unassociated with point contacts.  It can
        by used for joint stops however.
        """
        if name not in self._nslaws_data:
            nslaw=self._nslaws_data.create_dataset(name, (0,))
            nslaw.attrs['type']='RelayNSL'
            nslaw.attrs['size']=size
            nslaw.attrs['lb']=lb
            nslaw.attrs['ub']=ub
            nslaw.attrs['gid1']=collision_group1
            nslaw.attrs['gid2']=collision_group2

    def addJoint(self, name, object1, object2=None, pivot_point=[0, 0, 0],
                 axis=[0, 1, 0], joint_class='PivotJointR', absolute=None,
                 allow_self_collide=None, nslaws=None, stops=None, friction=None):
        """
        add a joint between two objects
        """
        if name not in self.joints():
            joint=self.joints().create_dataset(name, (0,))
            joint.attrs['object1']=object1
            if object2 is not None:
                joint.attrs['object2']=object2
            joint.attrs['type']=joint_class
            joint.attrs['pivot_point']=pivot_point
            joint.attrs['axis']=axis
            if absolute in [True, False]:
                joint.attrs['absolute']=absolute
            if allow_self_collide in [True, False]:
                joint.attrs['allow_self_collide']=allow_self_collide
            if nslaws is not None:
                joint.attrs['nslaws'] = nslaws # either name of one nslaw, or a
                                               # list of names same length as stops
            if stops is not None:
                joint.attrs['stops'] = stops # must be a table of [[axis,pos,dir]..]
            if friction is not None:
                joint.attrs['friction'] = friction # must be an NSL name (e.g.
                                                   # RelayNSL), or list of same


    def addBoundaryCondition(self, name, object1, indices=None, bc_class='HarmonicBC',
                             v=None, a=None, b=None, omega=None, phi=None):
        """
        add boundarycondition to the object object1

        implementation only works for HarmonicBC for the moment
        """
        if name not in self.boundary_conditions():
            boundary_condition=self.boundary_conditions().create_dataset(name, (0,))
            boundary_condition.attrs['object1']=object1
            boundary_condition.attrs['indices']=indices
            boundary_condition.attrs['type']=bc_class
            if bc_class == 'HarmonicBC' :
                boundary_condition.attrs['a']= a
                boundary_condition.attrs['b']= b
                boundary_condition.attrs['omega']= omega
                boundary_condition.attrs['phi']= phi
            elif bc_class == 'BoundaryCondition' :
                boundary_condition.attrs['v']= v
            elif bc_class == 'FixedBC' :
                pass # nothing to do
            else:
                raise NotImplementedError

    def run(self,
            with_timer=False,
            time_stepping=None,
            space_filter=None,
            options=None,
            body_class=None,
            shape_class=None,
            face_class=None,
            edge_class=None,
            controller=None,
            gravity_scale=1.0,
            t0=0,
            T=10,
            h=0.0005,
            multipoints_iterations=None,
            theta=0.50001,
            Newton_options=Kernel.SICONOS_TS_NONLINEAR,
            Newton_max_iter=20,
            set_external_forces=None,
            solver=Numerics.SICONOS_FRICTION_3D_NSGS,
            itermax=100000,
            tolerance=1e-8,
            projection_itermax=20,
            projection_tolerance=1e-8,
            projection_tolerance_unilateral=1e-8,
            numerics_verbose=False,
            violation_verbose=False,
            output_frequency=None,
            friction_contact_trace=False,
            friction_contact_trace_params=None,
            contact_index_set=1,
            osi=Kernel.MoreauJeanOSI):
        """
        Run a simulation from inputs in hdf5 file.
        parameters are:
          with_timer : use a timer for log output (default False)
          gravity_scale : gravity is multiplied by this factor.
                         1.     for meters (default).
                         1./100 for centimeters.
                         This parameter may be needed for small
                         objects because of Bullet collision margin (0.04).

          t0 : starting time (default 0)
          T  : end time      (default 10)
          h  : timestep      (default 0.0005)
          multiPointIterations : use bullet "multipoint iterations"
                                 (default True)
          theta : parameter for Moreau-Jean OSI (default 0.50001)
          Newton_max_iter : maximum number of iterations for
                          integrator Newton loop (default 20)
          set_external_forces : method for external forces
                                (default earth gravity)
          solver : default Numerics.SICONOS_FRICTION_3D_NSGS
          itermax : maximum number of iteration for solver
          tolerance : friction contact solver tolerance
          numerics_verbose : set verbose mode in numerics
          output_frequency :
          contact_index_set : index set from which contact points informations are retrieved.
        """
        print ('load siconos module ...')
        from siconos.kernel import \
            Model, NonSmoothDynamicalSystem, OneStepNSProblem,\
            TimeDiscretisation, GenericMechanical, FrictionContact,\
            NewtonImpactFrictionNSL

        from siconos.numerics import SICONOS_FRICTION_3D_ONECONTACT_NSN

        print ('setup model simulation ...')
        if set_external_forces is not None:
            self._set_external_forces=set_external_forces

        if space_filter is None:  space_filter = default_manager_class
        if time_stepping is None: time_stepping = default_simulation_class

        if output_frequency is not None:
            self._output_frequency=output_frequency

        if gravity_scale is not None:
            self._gravity_scale=gravity_scale

        # cold restart
        times=set()
        if self.dynamic_data() is not None and len(self.dynamic_data()) > 0:
            dpos_data=self.dynamic_data()
            times=set(dpos_data[:, 0])
            t0=float(max(times))

        # Time-related parameters for this simulation run
        k0=1+int(t0/h)
        k=k0
        kT=k0+int((T-t0)/h)
        if T > t0:
            print('')
            print('Simulation will run from {0:.4f} to {1:.4f}s, step {2} to step {3} (h={4}, times=[{5},{6}])'
                  .format(t0, T, k0, kT, h,
                          min(times) if len(times)>0 else '?',
                          max(times) if len(times)>0 else '?'))
            print('')
        else:
            print('Simulation time {0} >= T={1}, exiting.'.format(t0,T))
            exit()

        # Model
        #
        self._model=Model(t0, T)
        model=self._model

        # (1) OneStepIntegrators
        joints=list(self.joints())

        self._osi=osi(theta)

        # (2) Time discretisation --
        timedisc=TimeDiscretisation(t0, h)

        fc_index=0
        if (friction_contact_trace == False) :
            if len(joints) > 0:
                osnspb=GenericMechanical(SICONOS_FRICTION_3D_ONECONTACT_NSN)
                fc_index=1
            else:
                osnspb=FrictionContact(3, solver)
        else:
            osnspb=FrictionContactTrace(3, solver,friction_contact_trace_params,model)

        self._contact_index_set = contact_index_set

        # Global solver options
        solverOptions = osnspb.numericsSolverOptions()
        solverOptions.iparam[0]=itermax
        # -- full error evaluation
        #solverOptions.iparam[1]=Numerics.SICONOS_FRICTION_3D_NSGS_ERROR_EVALUATION_FULL
        # --  Adaptive error evaluation
        #solverOptions.iparam[1]=Numerics.SICONOS_FRICTION_3D_NSGS_ERROR_EVALUATION_ADAPTIVE
        #solverOptions.iparam[8]=1
        # -- light error evaluation with full final
        solverOptions.iparam[1] = Numerics.SICONOS_FRICTION_3D_NSGS_ERROR_EVALUATION_LIGHT
        solverOptions.iparam[14] = Numerics.SICONOS_FRICTION_3D_NSGS_FILTER_LOCAL_SOLUTION_TRUE
        solverOptions.dparam[0] = tolerance

        # Friction one-contact solver options
        fcOptions = solverOptions.internalSolvers[fc_index]
        fcOptions.solverId = Numerics.SICONOS_FRICTION_3D_ONECONTACT_NSN_GP_HYBRID
        fcOptions.iparam[0] = 100  # Local solver iterations

        osnspb.setMaxSize(30000)
        osnspb.setMStorageType(1)
        osnspb.setNumericsVerboseMode(numerics_verbose)

        # keep previous solution
        osnspb.setKeepLambdaAndYState(True)

        # Respect run() parameter for multipoints_iteratinos for
        # backwards compatibility, but this is overridden by
        # SiconosBulletOptions if one is provided.
        if multipoints_iterations is not None and options is None:
            options = SiconosBulletOptions()
            options.perturbationIterations = 3*multipoints_iterations
            options.minimumPointsPerturbationThreshold = 3*multipoints_iterations
        self._broadphase = space_filter(model, options)

        if use_original:
            if multipoints_iterations:
                if hasattr(self._broadphase, 'collisionConfiguration'):
                    self._broadphase.collisionConfiguration().\
                        setConvexConvexMultipointIterations()
                    self._broadphase.collisionConfiguration().\
                        setPlaneConvexMultipointIterations()
            else:
                print("""
            ConvexConvexMultipointIterations and PlaneConvexMultipointIterations are unset
            """)

        # (6) Simulation setup with (1) (2) (3) (4) (5)
        if time_stepping == Kernel.TimeSteppingDirectProjection:
            osnspb_pos=Kernel.MLCPProjectOnConstraints(Numerics.SICONOS_MLCP_ENUM, 1.0)
            so_pos = osnspb.numericsSolverOptions()
            so_pos.iparam[0]=itermax
            so_pos.dparam[0]=tolerance
            osnspb_pos.setMaxSize(30000)
            osnspb_pos.setMStorageType(0) # "not yet implemented for sparse storage"
            osnspb_pos.setNumericsVerboseMode(numerics_verbose)
            osnspb_pos.setKeepLambdaAndYState(True)
            simulation=time_stepping(timedisc, self._osi, osnspb, osnspb_pos)
            simulation.setProjectionMaxIteration(projection_itermax)
            simulation.setConstraintTolUnilateral(projection_tolerance_unilateral);
            simulation.setConstraintTol(projection_tolerance);
        else:
            simulation=time_stepping(timedisc)
            simulation.insertIntegrator(self._osi)
            simulation.insertNonSmoothProblem(osnspb)
        if use_proposed:
            simulation.insertInteractionManager(self._broadphase)

        simulation.setNewtonOptions(Newton_options)
        simulation.setNewtonMaxIteration(Newton_max_iter)
        simulation.setNewtonTolerance(1e-10)

        print ('import scene ...')
        self.importScene(t0, body_class, shape_class, face_class, edge_class)

        if controller is not None:
            controller.initialize(self)

        model.setSimulation(simulation)
        model.initialize()
        print ('first output static and dynamic objects ...')
        self.outputStaticObjects()
        self.outputDynamicObjects()

        if self._should_output_domains:
            log(self.outputDomains, with_timer)()

        # nsds=model.nonSmoothDynamicalSystem()
        # nds= nsds.getNumberOfDS()
        # for i in range(nds):
        #     ds=nsds.dynamicalSystem(i)
        #     ds.display()
        # raw_input()
        print ('start simulation ...')
        self._initializing=False
        while simulation.hasNextEvent():

            print ('step', k, '<', k0 + int((T - t0) / h))

            log(self.importBirths(body_class=body_class,
                                  shape_class=shape_class,
                                  face_class=face_class,
                                  edge_class=edge_class))

            if controller is not None:
                controller.step()

            if use_original:
                log(self._broadphase.buildInteractions, with_timer)\
                    (model.currentTime())

            if (friction_contact_trace == True) :
                 osnspb._stepcounter = k

            log(simulation.computeOneStep, with_timer)()

            if (k % self._output_frequency == 0) or (k == 1):
                print ('output in hdf5 file at step ', k)

                log(self.outputDynamicObjects, with_timer)()

                log(self.outputVelocities, with_timer)()

                log(self.outputContactForces, with_timer)()

                if self._should_output_domains:
                    log(self.outputDomains, with_timer)()

                log(self.outputSolverInfos, with_timer)()

                log(self._out.flush)()


            if use_proposed:
                number_of_contacts = (
                    self._broadphase.statistics().new_interactions_created
                    + self._broadphase.statistics().existing_interactions_processed)
            elif use_original:
                number_of_contacts = (self._model.simulation()
                                   .oneStepNSProblem(0).getSizeOutput()//3)

            if number_of_contacts > 0 :
                print('number of contacts',self._model.simulation().oneStepNSProblem(0).getSizeOutput()//3)
                self.printSolverInfos()

            if violation_verbose and number_of_contacts > 0 :
                if len(simulation.y(0,0)) >0 :
                    print('violation info')
                    y=simulation.y(0,0)
                    yplus=  np.zeros((2,len(y)))
                    yplus[0,:]=y
                    y=np.min(yplus,axis=1)
                    violation_max=np.max(-y)
                    print('  violation max :',violation_max)
                    if  (violation_max >= self._collision_margin):
                        print('  violation max is larger than the collision_margin')
                    lam=simulation.lambda_(1,0)
                    print('  lambda max :',np.max(lam))
                    #print(' lambda : ',lam)
                    #raw_input()


                if len(simulation.y(1,0)) >0 :
                    v=simulation.y(1,0)
                    vplus=  np.zeros((2,len(v)))
                    vplus[0,:]=v
                    v=np.max(vplus,axis=1)
                    print('  velocity max :',np.max(v))
                    print('  velocity min :',np.min(v))
                #     #print(simulation.output(1,0))


            log(simulation.nextStep, with_timer)()

            print ('')
            k += 1
