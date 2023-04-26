# Mechanics IO
"""Run a pre-configured Siconos "mechanics" HDF5 file."""

from __future__ import print_function

import os
import sys

from math import cos, sin, asin, atan2, acos
from scipy import constants

import numpy as np
import h5py
## fix compatibility with h5py version
if hasattr(h5py, 'vlen_dtype'):
    h5py_vlen_dtype = h5py.vlen_dtype
elif hasattr(h5py, 'new_vlen'):
    h5py_vlen_dtype = h5py.new_vlen

import bisect
import time
import numbers
import shutil
import json

import tempfile
from contextlib import contextmanager

# Siconos imports

import siconos.io.mechanics_hdf5
import siconos.numerics as sn
import siconos.kernel as sk
from siconos.kernel import \
    EqualityConditionNSL, \
    Interaction, DynamicalSystem, TimeStepping,\
    SICONOS_OSNSP_TS_VELOCITY,\
    cast_FrictionContact

# Siconos Mechanics imports
from siconos.mechanics.collision.tools import Contactor, Shape
from siconos.mechanics import joints
from siconos.io.io_base import MechanicsIO
from siconos.io.FrictionContactTrace import GlobalFrictionContactTrace as GFCTrace
from siconos.io.FrictionContactTrace import FrictionContactTrace as FCTrace
from siconos.io.FrictionContactTrace import GlobalRollingFrictionContactTrace as GRFCTrace
from siconos.io.mechanics_hdf5 import MechanicsHdf5


# Imports for mechanics 'collision' submodule
from siconos.mechanics.collision import RigidBodyDS, RigidBody2dDS, \
    SiconosSphere, SiconosBox, SiconosCylinder, SiconosCone, SiconosCapsule, \
    SiconosPlane, SiconosConvexHull, \
    SiconosDisk, SiconosBox2d, SiconosConvexHull2d, \
    SiconosContactor, SiconosContactorSet, \
    SiconosMesh, SiconosHeightMap

from siconos.mechanics.collision import SpaceFilter, bodies
from siconos.mechanics.collision.bodies import Disk,\
    Circle
from siconos.kernel import SimpleMatrix

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

# It is necessary to select a back-end, although currently only Bullet
# is supported for general objects.
backend = 'bullet'


def set_backend(b):
    global backend
    backend = b
    setup_default_classes()


have_bullet = False
have_occ = False

try:
    from siconos.mechanics.collision.bullet import \
        SiconosBulletCollisionManager, SiconosBulletOptions, SICONOS_BULLET_2D
    have_bullet = True
except Exception:
    have_bullet = False

# OCC imports
try:
    from siconos.mechanics import occ
    have_occ = True
except Exception:
    have_occ = False

# Configuration

default_manager_class = None
default_simulation_class = None
default_body_class = None
use_bullet = False


def setup_default_classes():
    global default_manager_class
    global default_simulation_class
    global default_body_class
    global use_bullet
    default_simulation_class = TimeStepping
    default_body_class = RigidBodyDS
    if backend == 'bullet':
        def m(bullet_options):
            if bullet_options is None:
                bullet_options = SiconosBulletOptions()
            return SiconosBulletCollisionManager(bullet_options)
        default_manager_class = m
        use_bullet = have_bullet
    elif backend == 'occ':
        default_manager_class = lambda options: occ.OccSpaceFilter()
        #def default_manager_class(options):
        #    occ.OccSpaceFilter()
        default_simulation_class = occ.OccTimeStepping
        default_body_class = occ.OccBody
        use_bullet = False
    elif backend == 'native':
        use_bullet = False
        def default_manager_class(options):
            sp =  SpaceFilter()
            sp.setBBoxfactor(3)
            sp.setCellsize(6)
            return sp


setup_default_classes()

# Constants

joint_points_axes = {
    'KneeJointR': (1, 0),
    'PivotJointR': (1, 1),
    'PrismaticJointR': (0, 1),
    'CylindricalJointR': (1, 1),
    'FixedJointR': (0, 0),
}


# Utility functions
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
def tmpfile(suffix='', prefix='siconos_io', contents=None,
            debug=False):
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
    if not debug:
        os.remove(tfilename)

time_measure = time.perf_counter
if (sys.version_info.major+0.1*sys.version_info.minor < 3.3):
    time_measure= time.clock



class Timer():

    def __init__(self):
        self._t0 = time_measure()

    def elapsed(self):
        return time_measure() - self._t0

    def update(self):
        self._t0 = time_measure()


def warn(msg):
    sys.stderr.write('{0}: {1}'.format(sys.argv[0], msg))


def object_id(obj):
    """returns an unique object identifier"""
    return obj.__hash__()


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
                                chunks=[None, (4000, nbcolumns)][comp],
                                compression=[None, 'gzip'][comp],
                                compression_opts=[None, 9][comp])

#
# misc fixes
#
# fix ctr.'name' in old hdf5 files
#
def upgrade_io_format(filename):

    with MechanicsHdf5(filename, mode='a') as io:

        for instance_name in io.instances():
            for contactor_instance_name in io.instances()[instance_name]:
                contactor = io.instances()[instance_name][
                    contactor_instance_name]
                if 'name' in contactor.attrs:
                    warn("""
contactor {0} attribute 'name': renamed in 'shape_name'
                    """)
                    contactor.attrs['shape_name'] = contactor['name']
                    del contactor['name']


def str_of_file(filename):
    with open(filename, 'r') as f:
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


#
# fix orientation -> rotation ?
#
def quaternion_get(orientation):
    """
    Get quaternion from orientation
    """
    if len(orientation) == 2:
        # axis + angle
        axis = orientation[0]
        if not (len(axis) == 3):
            raise AssertionError("quaternion_get. axis in not a 3D vector")
        angle = orientation[1]
        if not isinstance(angle, float):
            raise AssertionError("quaternion_get. angle must be a float")
        n = sin(angle / 2.) / np.linalg.norm(axis)
        ori = [cos(angle / 2.), axis[0] * n, axis[1] * n, axis[2] * n]
    else:
        if not (len(orientation) == 4):
            raise AssertionError(
                "quaternion_get. The quaternion must be of size 4")
        # a given quaternion
        ori = orientation
    return ori


def quaternion_multiply(q1, q0):
    w0, x0, y0, z0 = q0
    w1, x1, y1, z1 = q1
    return np.array([-x1 * x0 - y1 * y0 - z1 * z0 + w1 * w0,
                     x1 * w0 + y1 * z0 - z1 * y0 + w1 * x0,
                     -x1 * z0 + y1 * w0 + z1 * x0 + w1 * y0,
                     x1 * y0 - y1 * x0 + z1 * w0 + w1 * z0], dtype=np.float64)


def phi(q0, q1, q2, q3):
    """
    Euler angle phi from quaternion.
    """
    return atan2(2 * (q0 * q1 + q2 * q3), 1 - 2 * (q1 * q1 + q2 * q2))


def theta(q0, q1, q2, q3):
    """
    Euler angle theta from quaternion.
    """
    return asin(2 * (q0 * q2 - q3 * q1))


def psi(q0, q1, q2, q3):
    """
    Euler angle psi from quaternion.
    """
    return atan2(2 * (q0 * q3 + q1 * q2), 1 - 2 * (q2 * q2 + q3 * q3))


# vectorized versions
phiv = np.vectorize(phi)
thetav = np.vectorize(theta)
psiv = np.vectorize(psi)


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
        apoints = np.empty((3, num_points), dtype='f8')
        for i in range(0, points.GetNumberOfTuples()):
            p = points.GetTuple(i)
            apoints[0, i] = p[0]
            apoints[1, i] = p[1]
            apoints[2, i] = p[2]

        if scale is not None:
            apoints *= scale

        aindices = np.empty(num_triangles * 3, dtype=np.int)

        for i in range(0, num_triangles):
            c = polydata.GetCell(i)
            aindices[i * 3 + 0] = c.GetPointIds().GetId(0)
            aindices[i * 3 + 1] = c.GetPointIds().GetId(1)
            aindices[i * 3 + 2] = c.GetPointIds().GetId(2)

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

        if backend == 'native':
            self._primitive = {'Disk': NativeDiskShape,
                               'Circle': NativeCircleShape,
                               'Line': NativeLineShape}
        else:
            self._primitive = {'Sphere': SiconosSphere,
                               'Box': SiconosBox,
                               'Cylinder': SiconosCylinder,
                               'Capsule': SiconosCapsule,
                               'Cone': SiconosCone,
                               'Plane': SiconosPlane,
                               'Disk': SiconosDisk,
                               'Box2d': SiconosBox2d}

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

        if new_instance or shape_name not in self._shapes:

            # load shape if it is an existing file
            if not isinstance(self.url(shape_name), str) and \
               'primitive' not in self.attributes(shape_name):
                # assume a vtp file (xml) stored in a string buffer

                if self.attributes(shape_name)['type'] == 'vtp':
                    if self.shape(shape_name).dtype == h5py_vlen_dtype(str):
                        with tmpfile() as tmpf:
                            data = self.shape(shape_name)[:][0]
                            ## fix compatibility with h5py version: to be removed in the future
                            if (h5py.version.version_tuple.major >=3 ):
                                tmpf[0].write(data.decode("utf-8"))
                            else:
                                tmpf[0].write(data)
                            tmpf[0].flush()
                            scale = self.attributes(shape_name).get(
                                'scale', None)
                            mesh, dims = load_siconos_mesh(
                                tmpf[1], scale=scale)
                            self._shapes[shape_name] = mesh
                            mesh.setInsideMargin(
                                self.shape(shape_name).attrs.get(
                                    'insideMargin', min(dims) * 0.02))
                            mesh.setOutsideMargin(
                                self.shape(shape_name).attrs.get(
                                    'outsideMargin', 0))
                    else:
                        raise AssertionError(
                            'ShapeCollection.get(), attributes(shape_name)["type"] != vtp ')

                elif self.attributes(shape_name)['type'] in ['step', 'stp']:
                    from OCC.STEPControl import STEPControl_Reader
                    from OCC.BRep import BRep_Builder
                    from OCC.TopoDS import TopoDS_Compound
                    from OCC.IFSelect import IFSelect_RetDone,\
                        IFSelect_ItemsByEntity

                    builder = BRep_Builder()
                    comp = TopoDS_Compound()
                    builder.MakeCompound(comp)

                    if self.shape(shape_name).dtype != h5py_vlen_dtype(str) :
                        raise AssertionError("self.shape(shape_name).dtype != h5py_vlen_dtype(str)")

                    ## fix compatibility with h5py version: to be removed in the future
                    if (h5py.version.version_tuple.major >=3 ):
                        tmp_contents =self.shape(shape_name)[:][0].decode('utf-8')
                    else :
                        tmp_contents =self.shape(shape_name)[:][0]

                    with tmpfile(contents=tmp_contents) as tmpf:
                        step_reader = STEPControl_Reader()

                        status = step_reader.ReadFile(tmpf[1])

                        if status == IFSelect_RetDone:  # check status
                            failsonly = False
                            step_reader.PrintCheckLoad(
                                failsonly, IFSelect_ItemsByEntity)
                            step_reader.PrintCheckTransfer(
                                failsonly, IFSelect_ItemsByEntity)

                            #  ok = step_reader.TransferRoot(1)
                            step_reader.TransferRoots()
                            # VA : We decide to loads all shapes in the step file
                            nbs = step_reader.NbShapes()

                            for i in range(1, nbs + 1):
                                shape = step_reader.Shape(i)
                                builder.Add(comp, shape)

                            self._shapes[shape_name] = comp
                            self._io._keep.append(self._shapes[shape_name])

                elif self.attributes(shape_name)['type'] in ['iges', 'igs']:
                    from OCC.IGESControl import IGESControl_Reader
                    from OCC.BRep import BRep_Builder
                    from OCC.TopoDS import TopoDS_Compound
                    from OCC.IFSelect import IFSelect_RetDone,\
                        IFSelect_ItemsByEntity

                    builder = BRep_Builder()
                    comp = TopoDS_Compound()
                    builder.MakeCompound(comp)

                    if not (self.shape(shape_name).dtype == h5py_vlen_dtype(str)):
                        raise AssertionError("ShapeCollection.get()")

                    #assert(self.shape(shape_name).dtype == h5py_vlen_dtype(str))

                    with tmpfile(contents=self.shape(shape_name)[:][0]) as tmpf:
                        iges_reader = IGESControl_Reader()

                        status = iges_reader.ReadFile(tmpf[1])

                        if status == IFSelect_RetDone:  # check status
                            failsonly = False
                            iges_reader.PrintCheckLoad(
                                failsonly, IFSelect_ItemsByEntity)
                            iges_reader.PrintCheckTransfer(
                                failsonly, IFSelect_ItemsByEntity)

                            iges_reader.TransferRoots()
                            nbs = iges_reader.NbShapes()

                            for i in range(1, nbs + 1):
                                shape = iges_reader.Shape(i)
                                builder.Add(comp, shape)

                            self._shapes[shape_name] = comp
                            self._io._keep.append(self._shapes[shape_name])

                elif self.attributes(shape_name)['type'] in['brep']:
                    if 'contact' not in self.attributes(shape_name):

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
                        if not ('contact' in self.attributes(shape_name)):
                            raise AssertionError("contact not in self.attributes(shape_name)")
                        if not ('contact_index' in self.attributes(shape_name)):
                            raise AssertionError("contact_index not in self.attributes(shape_name)")
                        if not ('brep' in self.attributes(shape_name)):
                            raise AssertionError("brep not in self.attributes(shape_name)")


                        #assert 'contact' in self.attributes(shape_name)
                        #assert 'contact_index' in self.attributes(shape_name)
                        #assert 'brep' in self.attributes(shape_name)

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
                    # a heightmap defined by a matrix
                    hm_data = self.shape(shape_name)
                    r = hm_data.attrs['rect']
                    if (len(r) != 2):
                        raise AssertionError("len(r) != 2")
                    #assert(len(r) == 2)
                    hm = SiconosHeightMap(hm_data, r[0], r[1])
                    #dims = list(r) + [np.max(hm_data) - np.min(hm_data)]
                    #hm.setInsideMargin(
                    #    hm_data.attrs.get('insideMargin', np.min(dims) * 0.02))
                    hm.setInsideMargin(
                        hm_data.attrs.get('insideMargin', 0))
                    hm.setOutsideMargin(
                        hm_data.attrs.get('outsideMargin', 0))

                    self._shapes[shape_name] = hm

                elif self.attributes(shape_name)['type'] in ['convex']:
                    # a convex point set
                    points = self.shape(shape_name)
                    if self._io._dimension == 3:
                        convex = SiconosConvexHull(points)
                        dims = [points[:, 0].max() - points[:, 0].min(),
                                points[:, 1].max() - points[:, 1].min(),
                                points[:, 2].max() - points[:, 2].min()]
                    elif self._io._dimension == 2:
                        convex = SiconosConvexHull2d(points)
                        dims = [points[:, 0].max() - points[:, 0].min(),
                                points[:, 1].max() - points[:, 1].min()]

                    avoid_internal_edge_contact = self.shape(shape_name).attrs.get('avoid_internal_edge_contact', False)
                    if avoid_internal_edge_contact:
                        convex.setAvoidInternalEdgeContact(True)

                    convex.setInsideMargin(
                        self.shape(shape_name).attrs.get('insideMargin',
                                                         min(dims) * 0.02))
                    convex.setOutsideMargin(
                        self.shape(shape_name).attrs.get('outsideMargin', 0))
                    self._shapes[shape_name] = convex

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
                    box = primitive(attrs)
                    self._shapes[shape_name] = box
                    box.setInsideMargin(
                        self.shape(shape_name).attrs.get('insideMargin',
                                                         min(attrs)*0.02))
                    box.setOutsideMargin(
                        self.shape(shape_name).attrs.get('outsideMargin', 0))
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
                    prim.setInsideMargin(
                        shp.attrs.get('insideMargin', min(attrs)*0.02))
                    prim.setOutsideMargin(shp.attrs.get('outsideMargin', 0))

        return self._shapes[shape_name]


class MechanicsHdf5Runner_run_options(dict):
    def __init__(self):
        d={}
        d['with_timer']=False
        d['time_stepping']=None
        d['interaction_manager']=None
        d['bullet_options']=None
        d['body_class']=None
        d['shape_class']=None
        d['face_class']=None
        d['edge_class']=None
        d['controller']=None
        d['gravity_scale']=1.0
        d['t0']=0
        d['T']=10
        d['h']=0.0005
        d['multipoints_iterations']=None
        d['theta']=0.50001
        d['gamma']=0.0
        d['Newton_options']=sk.SICONOS_TS_NONLINEAR
        d['Newton_max_iter']=20
        d['Newton_tolerance']=1e-10
        d['set_external_forces']=None
        d['solver_options']=None
        d['solver_options_pos']=None
        d['osnspb_max_size']=0
        d['exit_tolerance']=None
        d['projection_itermax']=20
        d['projection_tolerance']=1e-8
        d['projection_tolerance_unilateral']=1e-8
        d['numerics_verbose']=False
        d['numerics_verbose_level']=0
        d['violation_verbose']=False
        d['verbose']=True
        d['verbose_progress']=True
        d['output_frequency']=None
        d['output_backup']=False
        d['output_backup_frequency']=None
        d['friction_contact_trace_params']=None
        d['output_contact_index_set']=1
        d['osi']=sk.MoreauJeanOSI
        d['constraint_activation_threshold']=0.0
        d['explode_Newton_solve']=False
        d['explode_computeOneStep']=False
        d['display_Newton_convergence']=False
        d['start_run_iteration_hook']=None
        d['end_run_iteration_hook']=None
        d['skip_last_update_output']=False
        d['skip_last_update_input']=False
        d['skip_reset_lambdas']=False
        d['osns_assembly_type']= None
        d['output_contact_forces']=True,
        d['output_contact_info']=True,


        super(self.__class__, self).__init__(d)


    def display(self):
        import pprint
        pp = pprint.PrettyPrinter(indent=4)
        pp.pprint(self)


class MechanicsHdf5Runner(siconos.io.mechanics_hdf5.MechanicsHdf5):

    """a Hdf5 context manager that reads the translations and
       orientations of collision objects from hdf5 file

       It provides functions to output translations and orientations in
       the same file during simulation (output is done by default in
       pos.dat)

       Parameters
       ----------
       io_filename: string, optional
            hdf5 file name, default = <caller>.hdf5, caller being the name
            without ext of the file that instanciates the Runner.
        mode: string, optional
            h5 mode (w, r, append), default = 'w'
        interaction_manager: SiconosCollisionManager, optional
            user-defined interaction handler (e.g. from Bullet), default=None
        nsds: NonSmoothDynamicalSystem, optional
            default = None
        simulation: Simulation, optional
            default = None
        osi: OneStepIntegrator, optional
            default = None
        shape_filename: string
            vtk file describing a mesh, default = None
        set_external_forces: python function, optional
            function used to apply forces onto the body.
            Must be :

            def funcname(body)

            body being a siconos Body (DS)
            Default : apply gravity forces.
        gravity_scale: real, optional
            multiplication factor for the gravity.
            1.     for meters (default).
            1./100 for centimeters.
            This parameter may be needed for small
            objects because of Bullet collision margin.
        collision_margin: real number, optional
            tolerance for collision, default = None (0.04 in Shape builder)
        use_compression: boolean, optional
            true to use compression for h5 file, default=False
        output_domains: boolean, optional
            if trueoutputs info regarding contact point domains
            default=False
        verbose: boolean, optional
           default=True

    """

    def __init__(self, io_filename=None, io_filename_backup=None, mode='w',
                 interaction_manager=None, nsds=None, simulation=None,
                 osi=None, shape_filename=None,
                 set_external_forces=None, gravity_scale=None,
                 collision_margin=None,
                 use_compression=False, output_domains=False, verbose=True):

        super(MechanicsHdf5Runner, self).__init__(io_filename, mode,
                                                  io_filename_backup,
                                                  use_compression,
                                                  output_domains, verbose)
        self._interman = interaction_manager
        self._nsds = nsds
        self._simulation = simulation
        self._osi = osi
        self._osnspb = None
        self._static = {}
        self._shape = None
        self._occ_contactors = dict()
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
        self._output_backup_frequency = 1
        self._keep = []
        self._scheduled_births = []
        self._scheduled_deaths = []
        self._births = dict()
        self._deaths = dict()
        self._initializing = True
        self._output_contact_index_set = 1
        self._start_run_iteration_hook = None
        self._end_run_iteration_hook = None
        self._ds_positions = None
        self._ds_boundary_conditions = {}
        self._time_stepping_class = None
        self._k = None
        self._k0 = None

    def __enter__(self):
        super(MechanicsHdf5Runner, self).__enter__()

        if self._gravity_scale is None:
            self._gravity_scale = 1  # 1 => m, 1/100. => cm

        if self._set_external_forces is None:
            self._set_external_forces = self.apply_gravity

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

    def log(self, fun, with_timer=False, before=True):
        if with_timer:
            t = Timer()

            def logged(*args):
                t.update()
                if before:
                    print('[io.mechanics]| {0:50s} ...'.format(fun.__name__),
                          end='', flush=True)
                else:
                    print('[io.mechanics]|-->start {0:42s} ...'.
                          format(fun.__name__), flush=True)
                output = fun(*args)
                endt = t.elapsed()
                if not before:
                    print('[io.mechanics]|-->end {0:44s} ...'.
                          format(fun.__name__), end='', flush=True)
                print('..... {0:6.2e} s'.format(endt))
                siconos.io.mechanics_hdf5.group(self.log_data(), fun.__name__)
                siconos.io.mechanics_hdf5.add_line(
                    siconos.io.mechanics_hdf5.data(self.log_data()
                                                   [fun.__name__],
                                                   'timing', 1), endt)
                if (isinstance(output, numbers.Number)):
                    siconos.io.mechanics_hdf5.add_line(
                        siconos.io.mechanics_hdf5.data(self.log_data()
                                                       [fun.__name__],
                                                       'value', 1),
                        float(output))
                return output
            return logged
        else:
            def silent(*args):
                output = fun(*args)
                return output
            return silent

    def apply_gravity(self, body):
        g = constants.g / self._gravity_scale
        if hasattr(body, 'scalarMass'):
            scalar_mass = body.scalarMass()
        elif hasattr(body, 'getMassValue'):
            scalar_mass = body.getMassValue()
        if self._dimension == 3:
            weight = [0, 0, - scalar_mass * g]
        elif self._dimension == 2:
            weight = [0, - scalar_mass * g, 0.]
        body.setFExtPtr(weight)

    def import_nonsmooth_law(self, name):
        if self._interman is not None:
            nslawClass = getattr(sk, self._nslaws_data[name].attrs['type'])
            if nslawClass == sk.NewtonImpactFrictionNSL:
                nslaw = nslawClass(
                    float(self._nslaws_data[name].attrs['e']), 0.,
                    float(self._nslaws_data[name].attrs['mu']),
                    self._dimension)
            elif nslawClass == sk.NewtonImpactRollingFrictionNSL:
                if self._dimension == 3:
                    nslaw = nslawClass(
                        float(self._nslaws_data[name].attrs['e']), 0.,
                        float(self._nslaws_data[name].attrs['mu']),
                        float(self._nslaws_data[name].attrs['mu_r']), 5)
                elif self._dimension == 2:
                    nslaw = nslawClass(
                        float(self._nslaws_data[name].attrs['e']), 0.,
                        float(self._nslaws_data[name].attrs['mu']),
                        float(self._nslaws_data[name].attrs['mu_r']), 3)
            elif nslawClass == sk.NewtonImpactNSL:
                nslaw = nslawClass(float(self._nslaws_data[name].attrs['e']))
            elif nslawClass == sk.RelayNSL:
                nslaw = nslawClass(int(self._nslaws_data[name].attrs['size']),
                                   float(self._nslaws_data[name].attrs['lb']),
                                   float(self._nslaws_data[name].attrs['ub']))
            if not nslaw:
                raise AssertionError("no nslaw")
            #assert(nslaw)
            self._nslaws[name] = nslaw
            gid1 = int(self._nslaws_data[name].attrs['gid1'])
            gid2 = int(self._nslaws_data[name].attrs['gid2'])
            if gid1 >= 0 and gid2 >= 0:
                self._interman.insertNonSmoothLaw(nslaw, gid1, gid2)


    def import_native_object(self, name, translation, orientation,
                             velocity, contactors, mass, given_inertia,
                             body_class, shape_class,
                             birth=False, number=None):

        if mass is None:
            # a static object
            # for native only plans
            self._static[name] = {
                'number': number,
                'origin': translation,
                'orientation': orientation}

            # only one contactor
            ctor = contactors[0]
            shape = self._shape.get(ctor.shape_name, new_instance=True)
            a = self._shape._io.shapes()[ctor.shape_name][:][0][0]
            b = self._shape._io.shapes()[ctor.shape_name][:][0][1]
            c = self._shape._io.shapes()[ctor.shape_name][:][0][2]

            self._interman.insertLine(a, b, c)

            body = None
            flag = 'static'
        else:
            # a dynamic object
            ctor = contactors[0]
            shape = self._shape.get(ctor.shape_name)
            attrs = self._shape.attributes(ctor.shape_name)

            if self._shape.attributes(ctor.shape_name)['primitive'] == 'Disk':
                body_class = Disk
            elif self._shape.attributes(ctor.shape_name)['primitive'] == 'Circle':
                body_class = Circle
            radius = self._shape._io.shapes()[ctor.shape_name][:][0][0]

            body = body_class(radius, mass,
                              list(translation)+list(orientation), velocity)

            self._set_external_forces(body)
            self._nsds.insertDynamicalSystem(body)
            if birth and self._verbose:
                self.print_verbose('birth of body named {0}, translation {1}, orientation {2}'.format(name, translation, orientation))
            flag = 'dynamic'

            if number is not None:
                body.setNumber(int(number))

        return body, flag

    def import_occ_object(self, name, translation, orientation,
                          velocity, contactors, mass, given_inertia,
                          body_class, shape_class, face_class, edge_class,
                          birth=False, number=None):


        if mass is None:
            # a static object
            body = None

            self._static[name] = {
                'number': number,
                'origin': translation,
                'orientation': orientation}
            flag='static'
        else:

            if not np.isscalar(mass) or mass <= 0:
                self.print_verbose('Warning mass must be a positive scalar')
                raise RuntimeError('Warning mass must be a positive scalar')

            if body_class is None:
                body_class = occ.OccBody

            #assert (given_inertia is not None)
            if given_inertia is  None:
                raise AssertionError("given_inertia is  None")
            inertia = given_inertia
            if inertia is not None:
                if np.shape(inertia) == (3,):
                    inertia = np.diag(inertia)
                elif np.shape(inertia) != (3, 3):
                    self.print_verbose('Wrong shape of inertia')

            body = body_class(
                list(translation) + list(orientation), velocity, mass, inertia)

            if number is not None:
                body.setNumber(int(number))
            flag='dynamic'

        ref_shape = {ctor.instance_name: occ.OccContactShape(
            self._shape.get(ctor.shape_name,
                            shape_class, face_class,
                            edge_class, new_instance=True))
                     for ctor in contactors}

        ref_added = dict()
        for contactor in contactors:

            contact_shape = None
            reference_shape = ref_shape[contactor.instance_name]

            self._keep.append(reference_shape)

            if hasattr(contactor, 'contact_type'):

                if contactor.contact_type == 'Face':
                    contact_shape = occ.OccContactFace(
                        reference_shape, contactor.contact_index)

                elif contactor.contact_type == 'Edge':
                    contact_shape = occ.OccContactEdge(
                        reference_shape, contactor.contact_index)

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

            if reference_shape not in ref_added:
                if body is not None:
                    body.addShape(reference_shape.shape(),
                                  contactor.translation,
                                  contactor.orientation)
                else:
                    occ.occ_move(
                        reference_shape.shape(),
                        list(contactor.translation) + list(contactor.orientation))
                ref_added[reference_shape] = True

        if body is not None:
            self._set_external_forces(body)

            # add the dynamical system to the non smooth
            # dynamical system
            self._nsds.insertDynamicalSystem(body)
            self._nsds.setName(body, str(name))

            if birth and self._verbose:
                self.print_verbose('birth of body named {0}, translation {1}, orientation {2}'.format(name, translation, orientation))

        return body, flag

    def import_bullet_object(self, name, translation, orientation,
                             velocity, contactors, mass, inertia,
                             body_class, shape_class, birth=False,
                             number=None):

        if body_class is None:
            body_class = default_body_class

        if self._dimension == 2:
            body_class = RigidBody2dDS

        if self._interman is not None and 'input' in self._data:
            body = None
            # ---------------
            # a static object
            # ---------------
            if mass is None:

                cset = SiconosContactorSet()
                csetpos = (translation + orientation)
                for c in contactors:
                    shp = self._shape.get(c.shape_name)
                    pos = list(c.translation) + list(c.orientation)
                    cset.append(SiconosContactor(shp, pos, c.group))
                    self.print_verbose('              Adding shape %s to static contactor'%c.shape_name, 'at relative position', pos)

                staticBody = self._interman.addStaticBody(cset, csetpos, number)
                self._static[name] = {
                    'number': number,
                    'origin': translation,
                    'orientation': orientation,
                    'shape': shp,
                }
                # In the case of a static object, we return the staticBody and a flag
                # this will be used to remove the contactor if we have a time of death
                return staticBody, 'static'

            # ---------------
            # a dynamic object
            # ---------------
            else:

                if not np.isscalar(mass) or mass <= 0:
                    self.print_verbose('Warning mass must be a positive scalar')

                inertia_ok = False
                if self._dimension == 3:
                    if inertia is not None:
                        if np.shape(inertia) == (3,):
                            inertia = np.diag(inertia)
                            inertia_ok = True
                        if np.shape(inertia) == (3, 3):
                            inertia_ok = True
                elif self._dimension == 2:
                    if inertia is not None:
                        if np.shape(inertia) == (1, 1) or\
                           np.shape(inertia) == (1, ) or np.isscalar(inertia):
                            inertia_ok = True

                if inertia_ok:
                    # create the dynamics object with mass and inertia
                    body = body_class(translation + orientation,
                                      velocity, mass, inertia)
                    body.setUseContactorInertia(False)
                else:
                    if inertia is not None:
                        if self._dimension == 3:
                            self.print_verbose('**** Warning inertia for object named {0} does not have the correct shape: {1} instead of (3, 3) or (3,)'.format(name, np.shape(inertia)))
                            self.print_verbose('**** Inertia will be computed with the shape of the first contactor')
                        elif self._dimension == 2:
                            self.print_verbose('**** Warning inertia for object named {0} does not have the correct shape: {1} instead of (1, 1) or (1,) or scalar'.format(name, np.shape(inertia)))
                            self.print_verbose('**** Inertia will be computed with the shape of the first contactor')

                    body = body_class(translation + orientation,
                                      velocity, mass)
                    body.setUseContactorInertia(True)


                fext = self._input[name].get('allow_self_collide',
                                             None)
                if fext is not None:
                    body.setFextPtr(fext)


                self_collide = self._input[name].get('allow_self_collide',
                                                     None)
                if self_collide is not None:
                    body.setAllowSelfCollide(not not self_collide)

                cset = SiconosContactorSet()
                for c in contactors:
                    shp = self._shape.get(c.shape_name)
                    pos = list(c.translation) + list(c.orientation)
                    cset.append(SiconosContactor(shp, pos, c.group))
                    self.print_verbose('              Adding shape %s to dynamic contactor'%c.shape_name, 'at relative position', pos)

                body.setContactors(cset)

            if body:
                # set id number
                if number is not None:
                    body.setNumber(int(number))

                # set external forces
                self._set_external_forces(body)

                # add the dynamical system to the non smooth
                # dynamical system
                self._nsds.insertDynamicalSystem(body)
                self._nsds.setName(body, str(name))

            return body, 'dynamic'

    def make_coupler_jointr(self, ds1_name, ds2_name, coupled, references):
        topo = self._nsds.topology()
        dof1, dof2, ratio = coupled[0, :]
        refds_name = None
        if len(references) > 0:
            joint1_name = references[0]
            joint1 = joints.cast_NewtonEulerJointR(
                topo.getInteraction(str(joint1_name)).relation())
            joint2 = joint1
            joint1_ds1 = self.joints()[joint1_name].attrs['object1']
            joint1_ds2 = self.joints()[joint1_name].attrs.get('object2', None)
        if len(references) > 1:
            joint2_name = references[1]
            joint2 = joints.cast_NewtonEulerJointR(
                topo.getInteraction(str(joint2_name)).relation())
            joint2_ds1 = self.joints()[joint2_name].attrs['object1']
            joint2_ds2 = self.joints()[joint2_name].attrs.get('object2', None)

            if len(references) > 2:
                refds_name = references[2]
            else:
                # if second joint provided but no reference, then
                # infer refds_name from the reference joints
                dss = set([joint1_ds1, joint1_ds2, joint2_ds1, joint2_ds2])
                diff = list(dss.difference(set([ds1_name, ds2_name])))
                # there must be exactly one reference in common that
                # is not either of the DSs
                refds_name = diff[0]

        if refds_name:
            refds = sk.cast_NewtonEulerDS(topo.getDynamicalSystem(str(refds_name)))

            # Determine reference indexes:
            # Assert if neither ds in reference joints is the
            # ref ds and that the other ds is ds1/ds2 as
            # appropriate.
            refidx1, refidx2 = 0, 0
            if joint1_ds1 == ds1_name:
                refidx1 = 1
                if (joint1_ds2 != refds_name):
                    raise AssertionError("joint1_ds2 != refds_name")
                #assert(joint1_ds2 == refds_name)
            else:
                if (joint1_ds1 != refds_name):
                    raise AssertionError("joint1_ds1 != refds_name")
                #assert(joint1_ds1 == refds_name)
                if (joint1_ds2 != ds1_name):
                    raise AssertionError("joint1_ds2 != ds1_name")
                #assert(joint1_ds2 == ds1_name)
            if joint2_ds1 == ds2_name:
                refidx2 = 1
                if (joint2_ds2 != refds_name):
                    raise AssertionError("joint2_ds2 != refds_name")
                #assert(joint2_ds2 == refds_name)
            else:
                #assert(joint2_ds1 == refds_name)
                if (joint2_ds1 != refds_name):
                    raise AssertionError("joint2_ds1 != refds_name")
                #assert(joint2_ds2 == ds2_name)
                if (joint2_ds2 != ds2_name):
                    raise AssertionError("joint2_ds2 != ds1_name")

            joint = joints.CouplerJointR(joint1, int(dof1),
                                         joint2, int(dof2), ratio,
                                         refds, refidx1, refds, refidx2)
        else:
            # otherwise we assume no reference, coupler is directly
            # between ds1 and ds2
            joint = joints.CouplerJointR(joint1, int(dof1),
                                         joint2, int(dof2), ratio)
        return joint

    def import_joint(self, name):
        if self._interman is not None:
            nsds = self._nsds
            topo = nsds.topology()

            joint_type = self.joints()[name].attrs['type']
            joint_class = getattr(joints, joint_type)
            absolute = not not self.joints()[name].attrs.get('absolute', True)
            allow_self_collide = self.joints()[name].attrs.get(
                'allow_self_collide', None)
            stops = self.joints()[name].attrs.get('stops', None)
            nslaws = self.joints()[name].attrs.get('nslaws', None)
            friction = self.joints()[name].attrs.get('friction', None)
            coupled = self.joints()[name].attrs.get('coupled', None)
            references = self.joints()[name].attrs.get('references', None)

            # work-around h5py unicode bug
            # https://github.com/h5py/h5py/issues/379
            if references is not None:
                references = [r.decode('utf-8') if isinstance(r, bytes) else r
                              for r in references]

            points = self.joints()[name].attrs.get('points', [])
            axes = self.joints()[name].attrs.get('axes', [])
            siconos.io.mechanics_hdf5.check_points_axes(name, joint_class,
                                                        points, axes)

            ds1_name = self.joints()[name].attrs['object1']
            ds1 = topo.getDynamicalSystem(ds1_name)
            ds2 = None

            if 'object2' in self.joints()[name].attrs:
                ds2_name = self.joints()[name].attrs['object2']
                ds2 = topo.getDynamicalSystem(ds2_name)

            if joint_class == joints.CouplerJointR:
                # This case is a little different, handle it specially
                #assert(references is not None)
                #assert(np.shape(coupled) == (1, 3))
                if (references is None):
                    raise AssertionError("references is None")
                if (np.shape(coupled) != (1, 3)):
                    raise AssertionError("np.shape(coupled) != (1, 3)")

                joint = self.make_coupler_jointr(ds1_name, ds2_name,
                                                 coupled, references)
                coupled = None  # Skip the use for "coupled" below, to
                                # install joint-local couplings

            else:
                # Generic NewtonEulerJointR interface
                joint = joint_class()
                for n, p in enumerate(points):
                    joint.setPoint(n, p)
                for n, a in enumerate(axes):
                    joint.setAxis(n, a)
                joint.setAbsolute(absolute)

            q1 = ds1.q()
            q2 = None if ds2 is None else ds2.q()

            joint.setBasePositions(q1, q2)

            if allow_self_collide is not None:
                joint.setAllowSelfCollide(not not allow_self_collide)
            joint_nslaw = EqualityConditionNSL(joint.numberOfConstraints())
            joint_inter = Interaction(joint_nslaw, joint)
            self._nsds.\
                link(joint_inter, ds1, ds2)
            nsds.setName(joint_inter, str(name))

            # Add a e=0 joint by default, otherwise user can specify
            # the impact law by name or a list of names for each axis.
            if stops is not None:
                if np.shape(stops)[1] != 3:
                    raise AssertionError("np.shape(stops)[1] != 3")
                #assert np.shape(stops)[1] == 3, 'Joint stops shape must be (?, 3)'
                if nslaws is None:
                    nslaws = [sk.NewtonImpactNSL(0.0)] * np.shape(stops)[0]
                elif isinstance(nslaws, bytes):
                    nslaws = [
                        self._nslaws[nslaws.decode('utf-8')]] * np.shape(stops)[0]
                elif isinstance(nslaws, str):
                    nslaws = [self._nslaws[nslaws]] * np.shape(stops)[0]
                else:
                    if np.shape(nslaws)[0] != np.shape(stops)[0]:
                        raise AssertionError("np.shape(nslaws)[0] != np.shape(stops)[0]")
                    #assert(np.shape(nslaws)[0] == np.shape(stops)[0])
                    nslaws = [self._nslaws[nsl] for nsl in nslaws]
                for n, (nsl, (axis, pos, _dir)) in enumerate(zip(nslaws, stops)):
                    # "bool()" is needed because type of _dir is
                    # numpy.bool_, which SWIG doesn't handle well.
                    stop = joints.JointStopR(joint, pos, bool(_dir < 0),
                                             int(axis))
                    stop_inter = Interaction(nsl, stop)
                    self._nsds.\
                        link(stop_inter, ds1, ds2)
                    nsds.setName(stop_inter, '%s_stop%d' % (str(name), n))

            # The per-axis friction NSL, can be ''
            if friction is not None:
                if isinstance(friction, str):
                    friction = [friction]
                elif isinstance(friction, bytes):
                    friction = [friction.decode('utf-8')]
                else:
                    friction = [(f.decode('utf-8') if isinstance(f, bytes) else f)
                                for f in friction]
                for ax, fr_nslaw in enumerate(friction):
                    if fr_nslaw == '':  # no None in hdf5, use empty string
                        continue        # instead for no NSL on an axis
                    nslaw = self._nslaws[fr_nslaw]
                    fr = joints.JointFrictionR(joint, [ax])
                    fr_inter = Interaction(nslaw, fr)
                    self._nsds.\
                        link(fr_inter, ds1, ds2)
                    nsds.setName(fr_inter, '%s_friction%d' % (str(name), ax))

            # An array of tuples (dof1, dof2, ratio) specifies
            # coupling between a joint's DoFs (e.g., to turn a
            # cylindrical joint into a screw joint)
            if coupled is not None:
                if len(coupled.shape) == 1:
                    coupled = np.array([coupled])
                if coupled.shape[1] != 3 :
                    raise AssertionError("coupled.shape[1] != 3")
                #assert(coupled.shape[1] == 3)
                for n, (dof1, dof2, ratio) in enumerate(coupled):
                    cpl = joints.CouplerJointR(joint, int(dof1),
                                               joint, int(dof2), ratio)
                    cpl.setBasePositions(q1, q2)
                    cpl_inter = Interaction(EqualityConditionNSL(1), cpl)
                    self._nsds.\
                        link(cpl_inter, ds1, ds2)
                    nsds.setName(cpl_inter, '%s_coupler%d' % (str(name), n))

    def import_boundary_conditions(self, name):
        if self._interman is not None:
            topo = self._nsds.\
                topology()

            bc_type = self.boundary_conditions()[name].attrs['type']
            bc_class = getattr(sk, bc_type)

            ds1_name = self.boundary_conditions()[name].attrs['object1']
            ds1 = topo.getDynamicalSystem(ds1_name)

            if bc_type == 'HarmonicBC':
                bc = bc_class(
                    self.boundary_conditions()[name].attrs['indices'],
                    self.boundary_conditions()[name].attrs['a'],
                    self.boundary_conditions()[name].attrs['b'],
                    self.boundary_conditions()[name].attrs['omega'],
                    self.boundary_conditions()[name].attrs['phi'])

            elif bc_type == 'FixedBC':
                bc = bc_class(
                    self.boundary_conditions()[name].attrs['indices'])

            elif(bc_type == 'BoundaryCondition'):
                bc = bc_class(
                    self.boundary_conditions()[name].attrs['indices'],
                    self.boundary_conditions()[name].attrs['v'])

            # set bc to the ds1

            ds1.setBoundaryConditions(bc)

            #joint_inter = Interaction(joint_nslaw, joint)
            #    self._nsds.\
            #        link(joint_inter, ds1)

    def import_permanent_interactions(self, name):
        """
        """
        if self._interman is not None and 'input' in self._data and\
           self.permanent_interactions() is not None:
            topo = self._nsds.topology()

            pinter = self.permanent_interactions()[name]
            body1_name = pinter.attrs['body1_name']
            body2_name = pinter.attrs['body2_name']

            try:
                ds1 = sk.cast_NewtonEulerDS(
                    topo.getDynamicalSystem(body1_name))
            except Exception:
                ds1 = None

            try:
                ds2 = sk.cast_NewtonEulerDS(
                    topo.getDynamicalSystem(body2_name))
            except Exception:
                ds2 = None

            # static object in second
            if ds1 is None:
                ds2, ds1 = ds1, ds2

            contactor1_name = pinter.attrs.get('contactor1_name', None)
            contactor2_name = pinter.attrs.get('contactor2_name', None)

            if contactor1_name is None:
                contactor1_names = self._occ_contactors[body1_name].keys()
            else:
                contactor1_names = [contactor1_name]

            if contactor2_name is None:
                contactor2_names = self._occ_contactors[body2_name].keys()
            else:
                contactor2_names = [contactor2_name]

            distance_calculator = pinter.attrs['distance_calculator']
            offset1 = pinter.attrs['offset1']
            offset2 = pinter.attrs['offset2']

            body1 = self._input[body1_name]
            body2 = self._input[body2_name]

            real_dist_calc = {'cadmbtb': occ.CadmbtbDistanceType,
                              'occ': occ.OccDistanceType}

            for contactor1_name in contactor1_names:
                for contactor2_name in contactor2_names:
                    ctr1 = body1[contactor1_name]
                    ctr2 = body2[contactor2_name]

                    cg1 = int(ctr1.attrs['group'])
                    cg2 = int(ctr2.attrs['group'])
                    nslaw = self._interman.nonSmoothLaw(cg1, cg2)

                    cocs1 = self._occ_contactors[body1_name][contactor1_name]
                    cocs2 = self._occ_contactors[body2_name][contactor2_name]

                    if ds2 is None:
                        self.print_verbose('moving contactor {0} of static object {1} to {2}'.format(contactor2_name, body2_name, list(np.array(body2.attrs['translation']) +\
                                                                                                                                       np.array(ctr2.attrs['translation'])) +\
                                                                                                     list(quaternion_multiply(ctr2.attrs['orientation'],
                                                                                                                              body2.attrs['orientation']))))
                        occ.occ_move(cocs2.data(),
                                     list(np.array(body2.attrs['translation']) +\
                                        np.array(ctr2.attrs['translation'])) +\
                                     list(quaternion_multiply(ctr2.attrs['orientation'],
                                          body2.attrs['orientation'])))

                    cp1 = occ.ContactPoint(cocs1)
                    cp2 = occ.ContactPoint(cocs2)

                    relation = occ.OccR(cp1, cp2,
                                        real_dist_calc[distance_calculator]())

                    relation.setOffset1(offset1)
                    relation.setOffset2(offset2)

                    inter = Interaction(nslaw, relation)

                    if ds2 is not None:
                        self._nsds.link(inter, ds1, ds2)
                    else:
                        self._nsds.link(inter, ds1)

                    # keep pointers
                    self._keep.append([cocs1, cocs2, cp1,
                                       cp2, relation])

    def import_object(self, name, body_class=None, shape_class=None,
                      face_class=None, edge_class=None, birth=False,
                      translation=None, orientation=None, velocity=None):
        """Import an object by name,
        possibly overriding initial position and velocity.
        """
        obj = self._input[name]
        self.print_verbose ('Import object name:', name)
        self.print_verbose ('              number (id): {0} '.format(obj.attrs['id']))

        if translation is None:
            translation = obj.attrs['translation']
        if orientation is None:
            orientation = obj.attrs['orientation']
        if velocity is None:
            velocity = obj.attrs['velocity']

        # bodyframe center of mass
        center_of_mass = floatv(obj.attrs.get('center_of_mass', [0, 0, 0]))

        mass = obj.attrs.get('mass', None)
        inertia = obj.attrs.get('inertia', None)

        if mass is None:
            self.print_verbose ('              static object')
            self.print_verbose ('              position', list(translation) + list(orientation))
        else:
            self.print_verbose ('              dynamic object')
            self.print_verbose ('              position', list(translation) + list(orientation))
            self.print_verbose ('              velocity', velocity)

        # self.print_verbose('mass = ', mass)
        # self.print_verbose('inertia = ', inertia)

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
                        shape_name=ctr.attrs['shape_name'],
                        collision_group=ctr.attrs['group'].astype(int),
                        contact_type=ctr.attrs['type'],
                        contact_index=ctr.attrs['contact_index'].astype(int),
                        relative_translation=np.subtract(ctr.attrs['translation'].astype(float),
                                                         center_of_mass),
                        relative_orientation=ctr.attrs['orientation'].astype(float)))
            elif 'group' in ctr.attrs:
                # bullet contact
                if occ_type :
                    raise AssertionError("occ type is found")
                #assert not occ_type
                contactors.append(
                    Contactor(
                        instance_name=ctr.attrs['instance_name'],
                        shape_name=ctr.attrs['shape_name'],
                        collision_group=ctr.attrs['group'].astype(int),
                        relative_translation=np.subtract(ctr.attrs['translation'].astype(float), center_of_mass),
                        relative_orientation=ctr.attrs['orientation'].astype(float)))
            else:
                # occ shape
                occ_type = True
                # fix: not only contactors here
                contactors.append(
                    Shape(
                        instance_name=ctr.attrs['instance_name'],
                        shape_name=ctr.attrs['shape_name'],
                        relative_translation=np.subtract(ctr.attrs['translation'].astype(float), center_of_mass),
                        relative_orientation=ctr.attrs['orientation'].astype(float)))

        if occ_type:
            # Occ object
            body, flag = self.import_occ_object(
                name, floatv(translation), floatv(orientation),
                floatv(velocity), contactors, mass,
                inertia, body_class, shape_class, face_class,
                edge_class, birth=birth,
                number=self.instances()[name].attrs['id'])
        elif backend == 'native':
            body, flag = self.import_native_object(
                name, floatv(translation), floatv(orientation),
                floatv(velocity), contactors, mass,
                inertia, body_class, shape_class, birth=birth,
                number=self.instances()[name].attrs['id'])
        else:
            # Bullet object
            body, flag = self.import_bullet_object(
                name, floatv(translation), floatv(orientation),
                floatv(velocity), contactors, mass,
                inertia, body_class, shape_class, birth=birth,
                number=self.instances()[name].attrs['id'])

        # import boundary conditions
        bc_name = self._ds_boundary_conditions.get(name, None)
        if bc_name is not None:
            self.import_boundary_conditions(bc_name)

        # schedule its death immediately
        time_of_death = obj.attrs.get('time_of_death', None)

        if time_of_death is not None :
            bisect.insort_left(self._scheduled_deaths, time_of_death)
            if time_of_death in self._deaths:
                self._deaths[time_of_death].append((name, obj, body, flag))
            else:
                self._deaths[time_of_death] = [(name, obj, body, flag)]

    def import_scene(self, time, body_class, shape_class, face_class,
                     edge_class):
        """
        From the specification given in the hdf5 file with the help of
        add* functions, import into the NSDS:
          - the static objects
          - the dynamic objects
          - the joints
        and into the interaction_manager:
          - the nonsmooth laws
        that have a specified time of birth <= current time.
        """

        # Ensure we count up from zero for implicit DS numbering
        DynamicalSystem.resetCount(0)

        for _ in self._ref:
            self._number_of_shapes += 1

        # import dynamical systems
        if self._interman is not None and 'input' in self._data:


            # We map the boundary conditions with the object into  self._ds_boundary_conditions.
            # In that way, boundary conditions will be imported once the
            # object is created in import_object and not necessarily at the
            # beginning of the simulation.
            for name in self.boundary_conditions():
                ds1_name = self.boundary_conditions()[name].attrs['object1']
                self._ds_boundary_conditions[ds1_name] = name
                # only one boundary condition per object. list is not needed
                # ds_bc = self._ds_boundary_conditions.get(ds1_name, None)
                # if ds_bc is None:
                #     self._ds_boundary_conditions[ds1_name] = [name]
                # else:
                #     self._ds_boundary_conditions[ds1_name].append(name)

            # get pointers
            dpos_data = self.dynamic_data()
            velocities = self.velocities_data()

            # some dict to prefetch values in order to
            # speedup cold start in the case of many objects
            xdpos_data = dict()
            xvelocities = dict()

            if dpos_data is not None and len(dpos_data) > 0:

                max_time = max(dpos_data[:, 0])
                id_last = np.where(
                    abs(dpos_data[:, 0] - max_time) < 1e-9)[0]
                id_vlast = np.where(abs(velocities[:, 0] - max_time) < 1e-9)[0]

                # prefetch positions and velocities in dictionaries
                # this avoids calls to np.where for each object
                for ldpos_data in dpos_data[id_last, :]:
                    xdpos_data[ldpos_data[1]] = ldpos_data
                for lvelocities in velocities[id_vlast, :]:
                    xvelocities[lvelocities[1]] = lvelocities

            else:
                # should not be used
                max_time = None
                id_last = None
            self.print_verbose('import dynamical systems ...')
            for (name, obj) in sorted(self._input.items(),
                                      key=lambda x: x[0]):

                mass = obj.attrs.get('mass', None)
                time_of_birth = obj.attrs.get('time_of_birth', -1)
                time_of_death = obj.attrs.get('time_of_death', float('inf'))

                if time_of_birth >= time:
                    #
                    # in the future
                    #
                    bisect.insort_left(self._scheduled_births, time_of_birth)
                    if time_of_birth in self._births:
                        self._births[time_of_birth].append((name, obj))
                    else:
                        self._births[time_of_birth] = [(name, obj)]
                elif time_of_death <= time:
                    # object already dead do not import
                    self.print_verbose('object', name, 'already dead do not import')
                    #input()
                    #pass
                else:
                    #
                    # this is for now
                    #
                    # cold restart if output previously done

                    if mass is not None and dpos_data is not None and\
                       len(dpos_data) > 0:
                        xpos = xdpos_data[obj.attrs['id']]
                        translation = (xpos[2], xpos[3], xpos[4])
                        orientation = (xpos[5], xpos[6], xpos[7], xpos[8])
                        xvel = xvelocities[obj.attrs['id']]
                        velocity = (xvel[2], xvel[3], xvel[4],
                                    xvel[5], xvel[6], xvel[7])
                        if self._dimension ==2:
                            angle=2.0*acos(translation[2])
                            orientation= (angle,)
                            translation = (translation[0],translation[1])
                            velocity =  (velocity[0],velocity[1],velocity[5])

                    else:
                        # start from initial conditions
                        translation = obj.attrs['translation']
                        orientation = obj.attrs['orientation']
                        velocity = obj.attrs['velocity']

                    self.import_object(name=name, body_class=body_class,
                                       shape_class=shape_class,
                                       face_class=face_class,
                                       edge_class=edge_class,
                                       translation=translation,
                                       orientation=orientation,
                                       velocity=velocity,
                                       birth=False)

            # import nslaws
            # note: no time of birth for nslaws and joints
            self.print_verbose('import nslaws ...')
            for name in self._nslaws_data:
                self.import_nonsmooth_law(name)


            # As for the  boundary conditions, the joints should be imported
            # after the importation of the object, or the two objects. If one
            # of the object has a time of birth this will fail.

            self.print_verbose('import joints ...')
            for name in self.joints():
                self.import_joint(name)

            self.print_verbose('import permanent interactions ...')
            for name in self.permanent_interactions():
                self.import_permanent_interactions(name)

    def current_time(self):
        if self._initializing:
            return self._simulation.startingTime()
        else:
            return self._simulation.nextTime()

    def import_births(self, body_class=None, shape_class=None,
                      face_class=None, edge_class=None,):
        """
        Import new objects into the NSDS.
        """
        time = self.current_time()

        ind_time = bisect.bisect_right(self._scheduled_births, time)

        current_times_of_births = set(self._scheduled_births[:ind_time])
        self._scheduled_births = self._scheduled_births[ind_time:]

        for time_of_birth in current_times_of_births:
            for (name, _) in self._births[time_of_birth]:

                self.import_object(name, body_class, shape_class,
                                   face_class, edge_class, birth=True)


    def execute_deaths(self):
        """
        Remove objects from the NSDS
        """
        time = self.current_time()

        ind_time = bisect.bisect_right(self._scheduled_deaths, time)
        #print('ind_time', ind_time)
        #print('self._scheduled_deaths', self._scheduled_deaths )
        #print('self._deaths', self._deaths )

        current_times_of_deaths = set(self._scheduled_deaths[:ind_time])
        #print('current_times_of_deaths', current_times_of_deaths )
        self._scheduled_deaths = self._scheduled_deaths[ind_time:]
        #print(self._scheduled_deaths)
        for time_of_death in current_times_of_deaths:
            #print(self._deaths[time_of_death])
            for ( _, _, body, flag) in self._deaths[time_of_death]:
                if flag == 'static':
                    self._interman.removeStaticBody(body)
                elif  flag == 'dynamic':
                    self._interman.removeBody(body)
                    self._nsds.removeDynamicalSystem(body)
                else:
                    msg = 'execute_deaths : unknown object type'
                    msg += 'It should static or dynamic'
                    raise RuntimeError(msg)
        #input()

    def output_static_objects(self):
        """
        Outputs translations and orientations of static objects
        """
        time = self.current_time()
        p = 0

        # append new static position
        current_line = self._static_data.shape[0]
        self._static_data.resize(current_line+len(self._static), 0)


        p=current_line
        for static in self._static.values():
            translation = static['origin']
            rotation = static['orientation']
            if self._dimension == 3:
                self._static_data[p, :] = \
                    [time,
                     static['number'],
                     translation[0],
                     translation[1],
                     translation[2],
                     rotation[0],
                     rotation[1],
                     rotation[2],
                     rotation[3]]
            elif self._dimension == 2:
                # VA. change the position such that is corresponds to a 3D object
                self._static_data[p, :] = \
                    [time,
                     static['number'],
                     translation[0],
                     translation[1],
                     0.0,
                     cos(rotation[0] / 2.0),
                     0.0,
                     0.0,
                     sin(rotation[0] / 2.0)]

            p += 1

        #print('current_line , self._static_data', current_line, self._static_data)
        #input()



    def output_dynamic_objects(self, initial=False):
        """
        Outputs translations and orientations of dynamic objects.
        """

        current_line = self._dynamic_data.shape[0]

        time = self.current_time()

        positions = self._io.positions(self._nsds)
        self._ds_positions = positions
        if positions is not None:
            self._dynamic_data.resize(current_line + positions.shape[0], 0)

            times = np.empty((positions.shape[0], 1))
            times.fill(time)
            if self._dimension == 3:
                self._dynamic_data[current_line:, :] = np.concatenate(
                    (times, positions), axis=1)
            elif self._dimension == 2:
                # VA. change the position such that is corresponds to a 3D object
                new_positions = np.zeros((positions.shape[0], 8))

                new_positions[:, 0] = positions[:, 0] # ds number
                new_positions[:, 1] = positions[:, 1] # x position
                new_positions[:, 2] = positions[:, 2] # y position

                new_positions[:, 4] = np.cos(positions[:, 3] / 2.0)
                new_positions[:, 7] = np.sin(positions[:, 3] / 2.0)
                self._dynamic_data[current_line:, :] = np.concatenate(
                    (times, new_positions), axis=1)

    def output_velocities(self):
        """
        Output velocities of dynamic objects
        """
        current_line = self._velocities_data.shape[0]

        time = self.current_time()

        velocities = self._io.velocities(self._nsds)

        if velocities is not None:

            self._velocities_data.resize(current_line + velocities.shape[0], 0)

            times = np.empty((velocities.shape[0], 1))
            times.fill(time)
            if self._dimension == 3:
                self._velocities_data[current_line:, :] = np.concatenate(
                    (times, velocities), axis=1)
            elif self._dimension == 2:
                 # VA. change the position such that is corresponds to a 3D object
                new_velocities = np.zeros((velocities.shape[0], 7))

                new_velocities[:, 0] = velocities[:, 0] # ds number
                new_velocities[:, 1] = velocities[:, 1] # x velocity
                new_velocities[:, 2] = velocities[:, 2] # y velocity

                new_velocities[:, 6] = velocities[:, 3] # theta velocity
                self._velocities_data[current_line:, :] = np.concatenate(
                    (times, new_velocities), axis=1)

    def output_contact_forces(self):
        """
        Outputs contact forces
        _output_contact_index_set default value is 1.
        """
        if self._nsds.\
                topology().indexSetsSize() > 1:
            time = self.current_time()
            contact_points = self._io.contactPoints(self._nsds,
                                                    self._output_contact_index_set)
            if contact_points is not None:
                current_line = self._cf_data.shape[0]
                # Increase the number of lines in cf_data
                # (h5 dataset with chunks)
                self._cf_data.resize(current_line + contact_points.shape[0], 0)
                times = np.empty((contact_points.shape[0], 1))
                times.fill(time)

                if self._dimension == 3:
                    self._cf_data[current_line:, :] = \
                        np.concatenate((times,
                                        contact_points),
                                       axis=1)

                elif self._dimension == 2:

                    # VA. change the contact info such that is corresponds to a 3D object
                    new_contact_points = np.zeros(
                        (contact_points.shape[0], 25))

                    new_contact_points[:, 0] = contact_points[:, 0]  # mu
                    new_contact_points[:, 1] = contact_points[:, 1]  # posa
                    new_contact_points[:, 2] = contact_points[:, 2]
                    #new_contact_points[:, 3]
                    new_contact_points[:, 4] = contact_points[:, 3]  # posb
                    new_contact_points[:, 5] = contact_points[:, 4]
                    #new_contact_points[:, 6]
                    new_contact_points[:, 7] = contact_points[:, 5]  # nc
                    new_contact_points[:, 8] = contact_points[:, 6]
                    #new_contact_points[:, 9]
                    # cf
                    new_contact_points[:, 10] = contact_points[:, 7]
                    new_contact_points[:, 11] = contact_points[:, 8]
                    #new_contact_points[:, 12]
                    new_contact_points[:, 13] = contact_points[:, 9]  # gap
                    new_contact_points[:, 14] = contact_points[:, 10]
                    #new_contact_points[:, 15]
                    # relative velocity
                    new_contact_points[:, 16] = contact_points[:, 11]
                    new_contact_points[:, 17] = contact_points[:, 12]
                    #new_contact_points[:, 18]
                    # reaction impulse
                    new_contact_points[:, 19] = contact_points[:, 13]
                    new_contact_points[:, 20] = contact_points[:, 14]
                    #new_contact_points[:, 21]
                    # inter id
                    new_contact_points[:, 22] = contact_points[:, 15]
                    new_contact_points[:, 23] = contact_points[:, 16]  # ds 1
                    new_contact_points[:, 24] = contact_points[:, 17]  # ds 2
                    self._cf_data[current_line:, :] = np.concatenate(
                        (times, new_contact_points), axis=1)

                # return the number of contacts
                return len(contact_points)
            return 0
        return 0
    def output_contact_info(self):
        """
        Outputs contact forces
        _output_contact_index_set default value is 1.
        """
        if self._nsds.\
                topology().indexSetsSize() > 1:
            time = self.current_time()
            contact_info = self._io.contactInfo(self._nsds,
                                                    self._output_contact_index_set)
            if contact_info is not None:
                current_line = self._cf_info.shape[0]
                # Increase the number of lines in cf_data
                # (h5 dataset with chunks)
                self._cf_info.resize(current_line + contact_info.shape[0], 0)
                times = np.empty((contact_info.shape[0], 1))
                times.fill(time)

                self._cf_info[current_line:, :] = \
                    np.concatenate((times,
                                    contact_info),
                                   axis=1)
                # return the number of contacts
                return len(contact_info)
            return 0
        return 0

    def output_domains(self):
        """
        Outputs domains of contact points
        """
        if self._nsds.\
                topology().indexSetsSize() > 1:
            time = self.current_time()
            domains = self._io.domains(self._nsds)

            if domains is not None:

                current_line = self._domain_data.shape[0]
                self._domain_data.resize(current_line + domains.shape[0], 0)
                times = np.empty((domains.shape[0], 1))
                times.fill(time)

                self._domain_data[current_line:, :] = \
                    np.concatenate((times, domains), axis=1)

    def output_solver_infos(self):
        """
        Outputs solver #iterations & precision reached
        """

        time = self.current_time()
        so = self._simulation.oneStepNSProblem(0).numericsSolverOptions()

        current_line = self._solv_data.shape[0]
        self._solv_data.resize(current_line + 1, 0)
        iterations = so.iparam[sn.SICONOS_IPARAM_ITER_DONE]
        precision = so.dparam[sn.SICONOS_DPARAM_RESIDU]
        if so.solverId == sn.SICONOS_GENERIC_MECHANICAL_NSGS:
            local_precision = so.dparam[3]  # Check this !
        elif so.solverId == sn.SICONOS_FRICTION_3D_NSGS:
            local_precision = 0.
        else:
            local_precision = precision

        self._solv_data[current_line, :] = [time, iterations, precision,
                                            local_precision]


    def output_results(self,with_timer=False):

        self.log(self.output_static_objects, with_timer)()

        self.log(self.output_dynamic_objects, with_timer)()

        self.log(self.output_velocities, with_timer)()

        if self._output_contact_forces:
            self.log(self.output_contact_forces, with_timer)()

        if self._output_contact_info and backend == 'bullet':
            self.log(self.output_contact_info, with_timer)()
        else:
            self.print_verbose('[warning] output_contact_info is only available with bullet backend for the moment')
            self.print_verbose('          to remove this message set output_contact_info options to False')

        if self._should_output_domains:
            self.log(self.output_domains, with_timer)()

        self.log(self.output_solver_infos, with_timer)()

        self.log(self._out.flush)()


    def output_run_options(self):
        """
        Outputs run_options
        """
        d= self._run_options.copy()

        so = d['solver_options']
        if so :
            d['solver_options'] = {}
            d['solver_options']['solverId'] = so.solverId
            d['solver_options']['solver name'] = so.solverId
            d['solver_options']['iparam size'] = so.iSize
            for i in range(so.iSize):
                d['solver_options']['iparam['+ str(i) +']'] = int(so.iparam[i])
            for i in range(so.dSize):
                if so.dparam[i] <= 1e24:
                    d['solver_options']['dparam['+ str(i)+']'] = float(so.dparam[i])
            # d['solver_options']['numberOfInternalSolvers']=so.numberOfInternalSolvers        # fix it

        sop = d['solver_options_pos']
        if sop :
            d['solver_options_pos'] = {}
            d['solver_options_pos']['solverId'] = so.solverId
            d['solver_options_pos']['solver name'] = so.solverId
            d['solver_options_pos']['iparam size'] = so.iSize
            for i in range(so.iSize):
                d['solver_options_pos']['iparam['+ str(i) +']'] = int(so.iparam[i])
            for i in range(so.dSize):
                if so.dparam[i] <= 1e24:
                    d['solver_options_pos']['dparam['+ str(i)+']'] = float(so.dparam[i])
            # d['solver_options_pos']['numberOfInternalSolvers']=so.numberOfInternalSolvers        # fix it

        bo = d['bullet_options']
        if bo:
            l = ['clearOverlappingPairCache', 'contactBreakingThreshold', 'contactProcessingThreshold',
                 'dimension', 'enablePolyhedralContactClipping', 'enableSatConvex',
                 'minimumPointsPerturbationThreshold', 'perturbationIterations',
                 'useAxisSweep3', 'worldScale']        # fix it
            l = ['contactProcessingThreshold', 'dimension', 'enablePolyhedralContactClipping',
                 'enableSatConvex', 'minimumPointsPerturbationThreshold', 'perturbationIterations',
                 'useAxisSweep3', 'worldScale']
            d['bullet_options'] = {}
            for e in l:
#                print('getattr(bo, e)', getattr(bo, e))
                d['bullet_options'][e] = getattr(bo, e)

        # to fix the serialization of run_options, we should use pickle
        # which is able to serialize python object.

        d['friction_contact_trace_params']='not serialized'   # fix it
        d['osi'] = 'not serialized'        # fix it
        d['time_stepping'] = 'not serialized'        # fix it

        d['start_run_iteration_hook']='not serialized'        # fix it
        d['end_run_iteration_hook']='not serialized'          # fix it
        if d['set_external_forces'] is not None:
            try :
                d['set_external_forces']= d['set_external_forces'].__name__ + '(name serialized)'
            except :
                d['set_external_forces']= 'not serialized'

        if d['controller'] is not None:
            try :
                d['controller']= d['controller'].__name__ + '(name serialized)'
            except :
                d['controller']= 'not serialized'

        dict_json=json.dumps(d)
        self._run_options_data.attrs['options'] = dict_json

    def print_solver_infos(self):
        """
        Outputs solver #iterations & precision reached
        """
        time = self.current_time()
        so = self._simulation.oneStepNSProblem(0).numericsSolverOptions()
        iterations = so.iparam[sn.SICONOS_IPARAM_ITER_DONE]
        precision = so.dparam[sn.SICONOS_DPARAM_RESIDU]
        msg = 'Numerics solver info at time : {0:10.6f}'.format(time)
        msg += ' iterations = {0:8d}'.format(iterations)
        msg += ' precision = {0:5.3e}'.format(precision)
        self.print_verbose(msg)

    def import_plugins(self):
        """
        Plugins extraction and compilation.
        """
        import subprocess

        for name in self._plugins:

            plugin_src = self._plugins[name][:]
            plugin_fname = self._plugins[name].attrs['filename']

            if os.path.isfile(plugin_fname):
                if str_of_file(plugin_fname) != plugin_src:
                    warn('plugin {0}: file {1} differs from the one stored hdf5'.format(name, plugin_fname))
            else:
                file_of_str(plugin_fname, plugin_src)

        # build
        subprocess.check_call(['siconos', '--build-plugins'])

    def import_external_functions(self):
        topo = self._nsds.topology()

        for name in self._external_functions:

            ext_fun = self._external_functions[name]
            plugin_name = ext_fun.attrs['plugin_name']
            plugin_function_name = ext_fun.attrs['plugin_function_name']
            body_name = ext_fun.attrs['body_name']

            ds = sk.cast_NewtonEulerDS(topo.getDynamicalSystem(body_name))

            if 'function_name' in ext_fun.attrs:
                function_name = ext_fun.attrs['function_name']

                getattr(ds, function_name)(plugin_name,
                                           plugin_function_name)
            else:
                bc_indices = ext_fun.attrs['bc_indices']
                # a boundary condition
                bc = sk.BoundaryCondition(bc_indices)
                bc.setComputePrescribedVelocityFunction(plugin_name,
                                                        plugin_function_name)
                ds.setBoundaryConditions(bc)

    def explode_Newton_solve(self, with_timer):
        s = self._simulation
        #self.log(s.initialize, with_timer)()
        self.log(s.initializeOSIAssociations, with_timer)()
        self.log(s.initializeIndexSets, with_timer)()

        self.log(s.updateWorldFromDS, with_timer)()
        self.log(s.updateInteractions, with_timer)()
        self.log(s.initializeNSDSChangelog, with_timer)()
        self.log(s.updateOutput, with_timer)()
        self.log(s.initOSNS, with_timer)()

        self.log(s.firstInitialize, with_timer)()
        if not s.skipResetLambdas():
            self.log(s.resetLambdas, with_timer)()
        # Again the access to s._newtonTolerance generates a segfault due to director
        #newtonTolerance = s.newtonTolerance()
        newtonMaxIteration = s.newtonMaxIteration()

        # return _kernel.TimeStepping_newtonSolve(self, criterion, maxStep)
        # RuntimeError: accessing protected member newtonSolve
        #s.newtonSolve(newtonTolerance, newtonMaxIteration);

        newtonNbIterations = 0
        isNewtonConverge = False
        explode_computeOneStep = self._run_options.get('explode_computeOneStep')

        #self.log(s.initializeNewtonLoop, with_timer)()
        if s.newtonOptions() == sk.SICONOS_TS_NONLINEAR :
            self.log(s.computeInitialNewtonState, with_timer)()
            self.log(s.computeResidu, with_timer)()
            self.log(s.updateWorldFromDS, with_timer)()
            self.log(s.updateInteractions, with_timer)()
            self.log(s.initializeNSDSChangelog, with_timer)()
            self.log(s.updateOutput, with_timer)()
            self.log(s.initOSNS, with_timer)()
            self.log(s.updateInput, with_timer)()

        self.log(s.updateDSPlugins, with_timer)(s.nextTime())
        self.log(s.computeResidu, with_timer)()

        # missing computeResiduY
        # self.log(s.computeResiduY, with_timer)()



        if s.newtonOptions() == sk.SICONOS_TS_LINEAR or s.newtonOptions() == sk.SICONOS_TS_LINEAR_IMPLICIT:
            if s.newtonOptions() == sk.SICONOS_TS_LINEAR_IMPLICIT:
                self.log(s.prepareNewtonIteration, with_timer)()
            self.log(s.computeFreeState, with_timer)()
            info=0
            if s.numberOfOSNSProblems() > 0:
                if explode_computeOneStep:
                    fc = self._osnspb
                    #self.log(fc.updateInteractionBlocks, with_timer)()
                    self.log(fc.preCompute, with_timer, before=False)(s.nextTime())
                    self.log(fc.updateMu, with_timer)()
                    if self._run_options.get('osi') == sk.MoreauJeanOSI :
                        if fc.getSizeOutput() != 0  :
                            info = self.log(fc.solve, with_timer)()
                            self.log(fc.postCompute, with_timer)()
                    else:
                        info = self.log(fc.solve, with_timer)()
                        self.log(fc.postCompute, with_timer)()
                else:
                    info = self.log(s.computeOneStepNSProblem, with_timer)(SICONOS_OSNSP_TS_VELOCITY)
                self.log(s.DefaultCheckSolverOutput, with_timer)(info)
                if not s.skipLastUpdateInput():
                    self.log(s.updateInput, with_timer)()
                self.log(s.updateState, with_timer)()
                if not s.skipLastUpdateOutput():
                    self.log(s.updateOutput, with_timer)()

        else:
            while (not isNewtonConverge) and newtonNbIterations < newtonMaxIteration:
                #self.print_verbose('newtonNbIterations',newtonNbIterations)
                info = 0
                newtonNbIterations = newtonNbIterations + 1
                self.log(s.prepareNewtonIteration, with_timer)()
                self.log(s.computeFreeState, with_timer)()
                if s.numberOfOSNSProblems() > 0:
                    if explode_computeOneStep:
                        fc = self._osnspb
                        self.log(fc.preCompute, with_timer)(s.nextTime())
                        self.log(fc.updateMu, with_timer)()
                        if fc.getSizeOutput() != 0 :
                            info = self.log(fc.solve, with_timer)()
                            self.log(fc.postCompute, with_timer)()
                    else:
                        info = self.log(s.computeOneStepNSProblem, with_timer)(SICONOS_OSNSP_TS_VELOCITY)
                self.log(s.DefaultCheckSolverOutput, with_timer)(info)
                self.log(s.updateInput, with_timer)()
                self.log(s.updateState, with_timer)()
                if (not isNewtonConverge) and (newtonNbIterations < newtonMaxIteration):
                    self.log(s.updateOutput, with_timer)()
                isNewtonConverge = self.log(s.newtonCheckConvergence, with_timer)
                if s.displayNewtonConvergence():
                    s.displayNewtonConvergenceInTheLoop()

            if s.displayNewtonConvergence():
                s.displayNewtonConvergenceAtTheEnd(info, newtonMaxIteration)

    def run_initialize(self,
            run_options=None,
            with_timer=False,
            time_stepping=None,
            interaction_manager=None,
            bullet_options=None,
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
            gamma=0.0,
            Newton_options=sk.SICONOS_TS_NONLINEAR,
            Newton_max_iter=20,
            set_external_forces=None,
            solver_options=None,
            solver_options_pos=None,
            osnspb_max_size=0,
            exit_tolerance=None,
            projection_itermax=20,
            projection_tolerance=1e-8,
            projection_tolerance_unilateral=1e-8,
            numerics_verbose=False,
            numerics_verbose_level=0,
            violation_verbose=False,
            verbose=True,
            verbose_progress=True,
            output_frequency=None,
            output_backup=False,
            output_backup_frequency=None,
            output_contact_forces=True,
            output_contact_info=True,
            friction_contact_trace_params=None,
            output_contact_index_set=1,
            osi=sk.MoreauJeanOSI,
            constraint_activation_threshold=0.0,
            explode_Newton_solve=False,
            explode_computeOneStep=False,
            display_Newton_convergence=False,
            start_run_iteration_hook=None,
            end_run_iteration_hook=None,
            skip_last_update_output=False
            ):
        """Run a simulation from a set of parameters described in a hdf5 file.

        Parameters
        ----------

        with_timer: boolean, optional
            if true, use a timer for log output (default False)
        time_stepping: siconos.kernel.Simulation, optional
             simulation type, default = sk.TimeStepping
        interaction_manager: SiconosCollisionManager, optional
            user-defined interaction handler (e.g. from Bullet), default=None
            (depends on the backend, e.g. Bullet or OCC).
            Warning: overwrite the value
            provided during MechanicsHdf5Runner init.
        bullet_options: ?, optional
            set of options for the interaction manager
            (e.g. SiconosBulletOptions), default = None
        body_class: siconos.mechanics.RigidBodyDS and heirs, optional
            class used for body definition, default = RigidBodyDS
        shape_class: ?, optional
            class used for shape definition (e.g. occ.OccContactShape)
            default = None
        face_class: ?, optional
            class used for face definition (e.g. occ.OccContactFace)
            default = None, (occ only?)
        edge_class=None,
            class used for edge definition (e.g. occ.OccContactEdge)
            default = None, (occ only?)
        controller: user-defined class, optional
            To apply user defined functions onto the nsds members.
            default = None, see example in wheels.py.
            The user-defined class must have two methods:

            initialize(self, MechanicsHdf5Runner)
            step()
        gravity_scale : int, optional
            multiplication factor for the gravity.
            1.     for meters (default).
            1./100 for centimeters.
            This parameter may be needed for small
            objects because of Bullet collision margin (0.04).
        t0: real, optional
            starting time (default 0)
        T: real, optional
            end time (default 10)
        h: real, optional
            time step size (default 0.0005)
        multiPoint_iterations : boolean, optional
            if true (default) use bullet "multipoint iterations"
        theta : real, optional
            parameter for Moreau-Jean OSI (default 0.50001)
        gamma : real, optional
            parameter for Moreau-Jean OSI (default 0.)
        Newton_options: int, optional
            sk.TimeStepping options to control the Newton loop
            (default sk.SICONOS_TS_NONLINEAR)
        Newton_max_iter : int, optional
            maximum number of iterations allowed for Newton integrator
            (default = 20).
        set_external_forces: python function, optional
            function used to apply forces onto the body.
            Must be :

            def funcname(body)

            body being a siconos Body (DS)
            Default : apply gravity forces.

            Note FP : already set in __enter__. Overwrite?

        solver_options : numerics SolverOptions, optional
            OneStepNsProblem solver options set.
            if solver_option is None, we leave Siconos/kernel choosing the default option
            (see solvers documentation for details)

        solver_options_pos : numerics SolverOptions for the position projection, optional
            OneStepNsProblem solver options set.
            if solver_option is None, we leave Siconos/kernel choosing the default option
            (see solvers documentation for details)

        osnspb_max_size : int, optional
            estimation of the maximum number of dynamical systems taken
            into account.
            Useful for memory pre-allocations and optimisation.
            if equal to 0 (default), it will be set to
            simulation().nonSmoothDynamicalSystem().topology().numberOfConstraints()

        exit_tolerance : real, optional
           if not None, the simulation will stop if precision >= exit_tolerance
           (default None).

        projection_itermax: int, optional
           max number of iteration for projection
           (only for TimeSteppingDirectProjection)
           default = 20
        projection_tolerance: real, optional
           tolerance for the violation of the equality constraints
           at the  position level (only for TimeSteppingDirectProjection)
           default = 1e-8
        projection_tolerance_unilateral=1e-8,
           tolerance for the violation of the unilateral constraints
           at the  position level (only for TimeSteppingDirectProjection)
           default = 1e-8

        numerics_verbose: boolean, optional
            True to activate numerics verbosity (default=False),
        numerics_verbose_level: int, optional,
            Set verbose level in numerics, default=0
        violation_verbose: boolean, optional
            If true, display information regarding constraint violation
            (if any), default=false
        verbose_progress: boolean, optional
            true to print current step informations (default=True)
        output_frequency : int, optional
            log and screen outputs frequency (default = 1)
        output_backup: boolean, optional
            True to backup hdf5 file (default false)
        output_backup_frequency: int, optional
            hdf5 file backup frequency (default = 1)
        friction_contact_trace_params: siconos.io.FrictionContactTraceParams,
            optional
            Set this to activate the wrapping of the one-step NS problem into
            FrictionContactTrace object. More log, more trace.
            Default = None.
        output_contact_index_set: int, optional
          index of the index set from which contact
          point information is retrieved. Default = 1
        osi: sk.OneStepIntegrator, optional
            class type used to describe one-step integration,
            default = sk.MoreauJeanOSI
        constraint_activation_threshold: real, optional
            threshold under which constraint is assume to be
            active. Default = 0.0,
        explode_Newton_solve: boolean, optional
            True to add more log/trace during Newton loop. Default=False,
        start_run_iteration_hook: boolean, optional
            if true, launch logging process at the beginning of each time step.
            Default = False.
        end_run_iteration_hook: boolean, optional
            if true, launch logging process at the end of each time step.
            Default = False.
        """

        if run_options is None:
            self.print_verbose('\nyou should consider to use a run_options dictionnary')
            self.print_verbose('we create a run_options for you and fill it with given options')
            self.print_verbose('some new options may not be available without run_options')

            run_options=MechanicsHdf5Runner_run_options()
            run_options['with_timer']=with_timer
            run_options['time_stepping']=time_stepping
            run_options['interaction_manager']=interaction_manager
            run_options['bullet_options']=bullet_options
            run_options['body_class']=body_class
            run_options['shape_class']=shape_class
            run_options['face_class']=face_class
            run_options['edge_class']=edge_class
            run_options['controller']=controller
            run_options['gravity_scale']=gravity_scale
            run_options['t0']=t0
            run_options['T']=T
            run_options['h']=h
            run_options['multipoints_iterations']=multipoints_iterations
            run_options['theta']=theta
            run_options['gamma']=gamma
            run_options['Newton_options']=Newton_options
            run_options['Newton_max_iter']=Newton_max_iter
            run_options['interaction_manager']= None
            run_options['set_external_forces']=set_external_forces
            run_options['solver_options']=solver_options
            run_options['solver_options_pos']=solver_options_pos
            run_options['osnspb_max_size']=osnspb_max_size
            run_options['exit_tolerance']=exit_tolerance
            run_options['projection_itermax']=projection_itermax
            run_options['projection_tolerance']=projection_tolerance
            run_options['projection_tolerance_unilateral']=projection_tolerance_unilateral
            run_options['numerics_verbose']=numerics_verbose
            run_options['numerics_verbose_level']=numerics_verbose_level
            run_options['violation_verbose']=violation_verbose
            run_options['verbose']=verbose
            run_options['verbose_progress']=verbose_progress
            run_options['output_frequency']=output_frequency
            run_options['output_backup']=output_backup
            run_options['output_backup_frequency']=output_backup_frequency
            run_options['friction_contact_trace_params']=friction_contact_trace_params
            run_options['output_contact_index_set']=output_contact_index_set
            run_options['osi']=osi
            run_options['constraint_activation_threshold']=constraint_activation_threshold
            run_options['explode_Newton_solve']=explode_Newton_solve
            run_options['explode_computeOneStep']=False,
            run_options['display_Newton_convergence']=display_Newton_convergence
            run_options['start_run_iteration_hook']=start_run_iteration_hook
            run_options['end_run_iteration_hook']=end_run_iteration_hook
            run_options['skip_last_update_output']=skip_last_update_output
            run_options['output_contact_forces']=output_contact_forces
            run_options['output_contact_info']=output_contact_info



        self._run_options=run_options

        self.print_verbose('run with run_options ...')

        if run_options['verbose']:
            run_options.display()

        self.print_verbose('setup model simulation ...')
        if run_options['set_external_forces'] is not None:
            self._set_external_forces = run_options['set_external_forces']

        interaction_manager=run_options['interaction_manager']
        if  interaction_manager is None:
            interaction_manager = default_manager_class

        if run_options['time_stepping'] is None:
            self._time_stepping_class = default_simulation_class
        else:
            self._time_stepping_class = run_options['time_stepping']

        if run_options['output_frequency'] is not None:
            self._output_frequency = run_options['output_frequency']

        if run_options['output_backup_frequency'] is not None:
            self._output_backup_frequency = run_options['output_backup_frequency']

        if run_options['output_backup'] is not None:
            self._output_backup = run_options['output_backup']

        if run_options['output_contact_forces'] is not None:
            self._output_contact_forces = run_options['output_contact_forces']

        if run_options['output_contact_info'] is not None:
            self._output_contact_info = run_options['output_contact_info']

        if run_options['gravity_scale'] is not None:
            self._gravity_scale = run_options['gravity_scale']

        t0=run_options['t0']
        T=run_options['T']
        h=run_options['h']

        # cold restart
        times = set()
        if self.dynamic_data() is not None and len(self.dynamic_data()) > 0:
            dpos_data = self.dynamic_data()
            times = set(dpos_data[:, 0])
            t0 = float(max(times))

        # Time-related parameters for this simulation run
        self._k0 = 1 + int(t0 / h)
        self._k = self._k0
        kT = self._k0 + int((T - t0) / h)
        if T > t0:
            self.print_verbose('')
            msg = 'Simulation will run from time {0:.4f} '.format(t0)
            msg += 'to {0:.4f}s, '.format(T)
            msg += 'step {} to step {} (h={}, '.format(self._k0, kT, h)
            msg += 'times=[{},{}])'.format(
                min(times) if len(times) > 0 else '?',
                max(times) if len(times) > 0 else '?')
            self.print_verbose(msg)
            self.print_verbose('')
        else:
            msg = 'Simulation time {0} >= T={1}, exiting.'.format(t0, T)
            self.print_verbose(msg)
            return


        # Respect run() parameter for multipoints_iterations for
        # backwards compatibility, but this is overridden by
        # SiconosBulletOptions if one is provided.
        multipoints_iterations =run_options.get('multipoints_iterations')
        bullet_options = run_options.get('bullet_options')
        if (multipoints_iterations is not None ) and (bullet_options is not None):
            msg = '[io.mechanics] run(): one cannot give multipoints_iterations and bullet_options simultaneously. \n'
            msg += '                             multipoints_iterations will be marked as obsolete. use preferably bullet_options\n'
            msg += '                             with  bullet_options.perturbationIterations and bullet_options.minimumPointsPerturbationThreshold.'
            raise RuntimeError(msg)

        if bullet_options is None and have_bullet:
            bullet_options = SiconosBulletOptions()
            if multipoints_iterations :
                bullet_options.perturbationIterations = 3 * multipoints_iterations
                bullet_options.minimumPointsPerturbationThreshold = \
                    3 * multipoints_iterations

        # MB: this may be in conflict with 'dimension' attribute
        if bullet_options is not None and \
           bullet_options.dimension == SICONOS_BULLET_2D:
            self._dimension = 2
        else:
            if (self._out.attrs.get('dimension', None) is None):
                # this is a second place to set the default
                self._dimension = 3

        self._interman = interaction_manager(bullet_options)


        joints = list(self.joints())
        if hasattr(self._interman, 'useEqualityConstraints'):
            if len(joints) == 0:
                self._interman.useEqualityConstraints(False)
            else:
                self._interman.useEqualityConstraints(True)

        # (0) NonSmooth Dynamical Systems definition
        self._nsds = sk.NonSmoothDynamicalSystem(t0, T)
        nsds = self._nsds

        self.print_verbose('import scene ...')
        body_class=run_options.get('body_class')
        shape_class=run_options.get('shape_class')
        face_class=run_options.get('face_class')
        edge_class=run_options.get('egde_class')

        self.import_scene(t0, body_class, shape_class, face_class, edge_class)

        self._output_contact_index_set = run_options.get('output_contact_index_set')

        # (1) OneStepIntegrators
        osi = run_options.get('osi')
        self._osi = osi(run_options.get('theta'))
        self._osi.setConstraintActivationThreshold(
            run_options['constraint_activation_threshold'])

        if run_options.get('gamma'):
            self._osi.setGamma(run_options.get('gamma'))


        # (2) Time discretisation --
        timedisc = sk.TimeDiscretisation(t0, h)

        # (3) choice of default OneStepNonSmoothProblem
        # w.r.t the type of nslaws
        nslaw_type_list = []
        for name in self._nslaws_data:
            nslaw_type_list.append(self._nslaws_data[name].attrs['type'])

        #print(set(nslaw_type_list))

        # This trick is used to add the EqualityConditionNSL
        # to the list of nslaw type
        # this must be improved by adding the EqualityConditionNSL
        # in self._nslaws_data
        # when a joint is imported.
        # For the moment, the nslaw is implicitely added
        # when we import_joint but is not stored
        # self._nslaws_data

        if len(joints) > 0:
            nslaw_type_list.append('EqualityConditionNSL')

        nb_of_nslaw_type = len(set(nslaw_type_list))

        # Check nslaw/OSI compatibility
        if (osi == sk.MoreauJeanGOSI):
            checknsl = 'NewtonImpactFrictionNSL' in set(nslaw_type_list) or 'NewtonImpactRollingFrictionNSL' in set(nslaw_type_list)
            if  not checknsl:
                msg = 'MoreauJeanGOSI can only deal '
                msg += 'with NewtonImpactFrictionNSL or NewtonImpactRollingFrictionNSL.'
                raise RuntimeError(msg)
            if (nb_of_nslaw_type > 1):
                msg = 'MoreauJeanGOSI cannot only deal with multiple inpact laws at the same time '
                raise RuntimeError(msg)

        # Creates and initialises the one-step nonsmooth problem.
        # The OSI and/or the nonsmooth law drives the choice.

        # rationale for choosing numerics solver options
        # if solver_option is None --> we leave Siconos/kernel choosing the default option
        # else we use the user solver_options
        solver_options= run_options.get('solver_options')
        osnspb_max_size = run_options.get('osnspb_max_size')
        osnspb_assembly_type =  run_options.get('osns_assembly_type')
        friction_contact_trace_params = run_options.get('friction_contact_trace_params')
        if friction_contact_trace_params is None:
            # Global friction contact.
            if (osi == sk.MoreauJeanGOSI):
                if 'NewtonImpactFrictionNSL' in set(nslaw_type_list):
                    if (solver_options is None):
                        osnspb = sk.GlobalFrictionContact(self._dimension)
                    else:
                        osnspb = sk.GlobalFrictionContact(self._dimension, solver_options)
                elif 'NewtonImpactRollingFrictionNSL' in set(nslaw_type_list):
                    if self._dimension ==3:
                        dimension_contact=5
                    elif self._dimension ==2:
                        dimension_contact=3
                    if (solver_options is None):
                        osnspb = sk.GlobalRollingFrictionContact(dimension_contact)
                    else:
                        osnspb = sk.GlobalRollingFrictionContact(dimension_contact, solver_options)
                osnspb.setMStorageType(sn.NM_SPARSE)
                # if sid == sn.SICONOS_GLOBAL_FRICTION_3D_ADMM:
                #     osnspb.setMStorageType(sn.NM_SPARSE)
                #     # which is the default for gfc in kernel, is it?
                # else:
                #     osnspb.setMStorageType(sn.NM_SPARSE_BLOCK)
                osnspb.setMaxSize(osnspb_max_size)
            else:
                if ('EqualityConditionNSL' in set(nslaw_type_list)):
                    if (solver_options is None):
                        osnspb = sk.GenericMechanical()
                    else:
                        osnspb = sk.GenericMechanical(solver_options)
                else:
                    if ('NewtonImpactFrictionNSL' in set(nslaw_type_list)) or \
                       (len(set(nslaw_type_list)) == 0):
                        if (solver_options is None):
                            osnspb = sk.FrictionContact(self._dimension)
                        else:
                            osnspb = sk.FrictionContact(self._dimension,solver_options)
                        osnspb.setMaxSize(osnspb_max_size)
                        osnspb.setMStorageType(sn.NM_SPARSE_BLOCK)


                    elif 'NewtonImpactRollingFrictionNSL' in set(nslaw_type_list):
                        if self._dimension ==3:
                            dimension_contact=5
                        elif self._dimension ==2:
                            dimension_contact=3
                        if (solver_options is None):
                            osnspb = sk.RollingFrictionContact(dimension_contact)
                        else:
                            osnspb = sk.RollingFrictionContact(dimension_contact, solver_options)
                    else:
                        msg = "Unknown nslaw type"
                        msg += str(set(nslaw_type_list))
                        raise RuntimeError(msg)

            if osnspb_assembly_type :
                osnspb.setMStorageType(sn.NM_SPARSE)
                osnspb.setAssemblyType(osnspb_assembly_type)

        else:  # With trace
            if solver_options is None:
                solver_options = sk.solver_options_create(
                    sn.SICONOS_FRICTION_3D_NSGS)
            sid = solver_options.solverId
            if(osi == sk.MoreauJeanGOSI):
                if 'NewtonImpactFrictionNSL' in set(nslaw_type_list):
                    osnspb = GFCTrace(3, solver_options, friction_contact_trace_params, nsds)
                    osnspb.setMStorageType(sn.NM_SPARSE)
                    osnspb.setMaxSize(osnspb_max_size)
                elif 'NewtonImpactRollingFrictionNSL' in set(nslaw_type_list):
                    osnspb = GRFCTrace(5, solver_options, friction_contact_trace_params, nsds)
                    osnspb.setMStorageType(sn.NM_SPARSE)
                    osnspb.setMaxSize(osnspb_max_size)
                else:
                    msg = "Unknown nslaw type"
                    msg += str(set(nslaw_type_list))
                    raise RuntimeError(msg)
            else:
                osnspb = FCTrace(3, solver_options, friction_contact_trace_params, nsds)
                osnspb.setMaxSize(osnspb_max_size)
                osnspb.setMStorageType(sn.NM_SPARSE_BLOCK)

        numerics_verbose=run_options.get('numerics_verbose')
        osnspb.setNumericsVerboseMode(numerics_verbose)
        if numerics_verbose:
            sn.numerics_set_verbose(run_options.get('numerics_verbose_level'))

        # keep previous solution
        osnspb.setKeepLambdaAndYState(True)


        self._osnspb = osnspb


        # (6) Simulation setup with (1) (2) (3) (4) (5)
        if self._time_stepping_class == sk.TimeSteppingDirectProjection:
            if solver_options_pos is None:
                osnspb_pos = sk.MLCPProjectOnConstraints(sn.SICONOS_MLCP_ENUM, 1.0)
            else:
                osnspb_pos = sk.MLCPProjectOnConstraints(solver_options_pos, 1.0)

            osnspb_pos.setMaxSize(osnspb_max_size)
            osnspb_pos.setMStorageType(sn.NM_DENSE)
            # "not yet implemented for sparse storage"
            osnspb_pos.setNumericsVerboseMode(numerics_verbose)
            osnspb_pos.setKeepLambdaAndYState(True)
            simulation = self._time_stepping_class(nsds, timedisc, self._osi,
                                       osnspb, osnspb_pos)
            simulation.setProjectionMaxIteration(projection_itermax)
            simulation.setConstraintTolUnilateral(
                projection_tolerance_unilateral)
            simulation.setConstraintTol(projection_tolerance)
        else:
            simulation = self._time_stepping_class(nsds, timedisc)
            simulation.insertIntegrator(self._osi)
            simulation.insertNonSmoothProblem(osnspb)

        simulation.insertInteractionManager(self._interman)

        simulation.setNewtonOptions(run_options['Newton_options'])
        simulation.setNewtonMaxIteration(run_options['Newton_max_iter'])
        simulation.setNewtonTolerance(run_options['Newton_tolerance'])


        simulation.setSkipLastUpdateOutput(run_options.get('skip_last_update_output'))
        simulation.setSkipLastUpdateInput(run_options.get('skip_last_update_input'))
        simulation.setSkipResetLambdas(run_options.get('skip_reset_lambdas'))

        verbose=run_options.get('verbose')
        if verbose:
            simulation.setDisplayNewtonConvergence(run_options.get('display_Newton_convergence'))

        self._simulation = simulation

        if len(self._plugins) > 0:
            self.print_verbose('import plugins ...')
            self.import_plugins()

        if len(self._external_functions) > 0:
            self.print_verbose('import external functions ...')
            self.import_external_functions()
        controller = run_options.get('controller')
        if controller is not None:
            controller.initialize(self)

        self._start_run_iteration_hook= run_options.get('start_run_iteration_hook')
        if self._start_run_iteration_hook is not None:
            self._start_run_iteration_hook.initialize(self)
        self._end_run_iteration_hook= run_options.get('end_run_iteration_hook')
        if self._end_run_iteration_hook is not None:
            self._end_run_iteration_hook.initialize(self)

        self.print_verbose('first output static and dynamic objects ...')
        self.output_static_objects()
        self.output_dynamic_objects()
        self.output_velocities()


        if self._should_output_domains:
            self.log(self.output_domains, with_timer)()


        self.output_run_options()

        # nsds = model.nonSmoothDynamicalSystem()
        # nds= nsds.getNumberOfDS()
        # for i in range(nds):
        #     ds = nsds.dynamicalSystem(i)
        #     ds.display()

        self.print_verbose('start simulation ...')
        self._initializing = False

    def run_loop(self):
        verbose = self._run_options.get('verbose')
        self._verbose=verbose
        with_timer= self._run_options.get('with_timer')
        body_class = self._run_options.get('body_class')
        shape_class = self._run_options.get('shape_class')
        face_class = self._run_options.get('face_class')
        edge_class = self._run_options.get('edge_class')
        controller = self._run_options.get('controller')
        friction_contact_trace_params = self._run_options.get('friction_contact_trace_params')
        exit_tolerance = self._run_options.get('exit_tolerance')
        t0 = self._run_options.get('t0')
        T = self._run_options.get('T')
        h = self._run_options.get('h')
        while self._simulation.hasNextEvent():

            if self._run_options.get('verbose_progress'):
                self.print_verbose('step', self._k, 'of', self._k0 + int((T - t0) / h) - 1)

            if self._start_run_iteration_hook is not None:
                if False == self.log(self._start_run_iteration_hook.call, with_timer, before=False)(self._k):
                    break

            self.log(self.import_births, with_timer)(body_class,
                                                     shape_class,
                                                     face_class,
                                                     edge_class)

            self.log(self.execute_deaths, with_timer)()

            if controller is not None:
                controller.step()

            if friction_contact_trace_params is not None:
                osnspb._stepcounter = self._k

            if self._run_options.get('explode_Newton_solve'):
                if(self._time_stepping_class == sk.TimeStepping):
                    self.log(self.explode_Newton_solve, with_timer,
                             before=False)(with_timer)
                else:
                    self.print_verbose('| [warning]. simulation of type', self._time_stepping_class.__name__,
                          ' has no exploded version')
                    self.log(self._simulation.computeOneStep, with_timer)()
            else:
                self.log(self._simulation.computeOneStep, with_timer)()

            cond = self._output_frequency and (self._k % self._output_frequency == 0)
            if cond or self._k == 1:
                if verbose:
                    self.print_verbose('output in hdf5 file at step ', self._k, ' time =', self.current_time())

                self.log(self.output_results, with_timer)()

            if self._output_backup:
                if (self._k % self._output_backup_frequency == 0) or (self._k == 1):

                    # close io file, hdf5 memory is cleaned
                    self._out.close()
                    try:
                        shutil.copyfile(self._io_filename,
                                        self._io_filename_backup)
                    except shutil.Error as e:
                        warn(str(e))
                    # open the file again
                    finally:
                        self.__enter__()

            self.log(self._simulation.clearNSDSChangeLog, with_timer)()

            # Note these are not the same and neither is correct.
            # "_interman.statistics" gives the number of contacts
            # collected by the collision engine, but it's possible some
            # are not in indexset1.  Meanwhile checking the size of
            # the non-smooth problem is wrong when there are joints.
            if use_bullet:
                number_of_contacts = self._interman.statistics().new_interactions_created
                number_of_contacts += self._interman.statistics().existing_interactions_processed
                if verbose and number_of_contacts > 0:
                    bullet_statistics = self._interman.statistics()
                    self.print_verbose('bullet_statistics:',
                                       'new_interactions_created :', bullet_statistics.new_interactions_created,
                                       'existing_interactions_processed :', bullet_statistics.existing_interactions_processed,
                                       'interaction_warnings :', bullet_statistics.interaction_warnings)
                    self.print_verbose('number of contacts',
                                       number_of_contacts,
                                       '(detected)',
                                       self._osnspb.getSizeOutput() // self._dimension ,
                                       '(active at velocity level. approx)')
                    self.print_solver_infos()

            else:
                number_of_contacts = self._osnspb.getSizeOutput() // self._dimension
                if verbose and number_of_contacts > 0:
                    msg = 'number of active contacts at the velocity level (approx)'
                    self.print_verbose(msg, number_of_contacts)
                    self.print_solver_infos()

            if self._run_options.get('violation_verbose') and number_of_contacts > 0:
                if len(self._simulation.y(0, 0)) > 0:
                    self.print_verbose('violation info')
                    y = self._simulation.y(0, 0)
                    yplus = np.zeros((2, len(y)))
                    yplus[0, :] = y
                    y = np.min(yplus, axis=1)
                    violation_max = np.max(-y)
                    self.print_verbose('  violation max :', violation_max)
                    if self._collision_margin is not None:
                        if(violation_max >= self._collision_margin):
                            self.print_verbose('  violation max is larger than the collision_margin')
                    lam = self._simulation.lambda_(1, 0)
                    self.print_verbose('  lambda max :', np.max(lam))
                    #print(' lambda : ',lam)

                if len(self._simulation.y(1, 0)) > 0:
                    v = self._simulation.y(1, 0)
                    vplus = np.zeros((2, len(v)))
                    vplus[0, :] = v
                    v = np.max(vplus, axis=1)
                    self.print_verbose('  velocity max :', np.max(v))
                    self.print_verbose('  velocity min :', np.min(v))
                #     #print(self._simulation.output(1,0))

            if (exit_tolerance is not None):
                solver_options = self._osnspb.numericsSolverOptions()
                precision = solver_options.dparam[sn.SICONOS_DPARAM_RESIDU]
                if (precision > exit_tolerance):
                    print('precision is larger exit_tolerance')
                    return False

            self.log(self._simulation.nextStep, with_timer)()

            if self._end_run_iteration_hook is not None:
                if False == self.log(self._end_run_iteration_hook.call, with_timer)(self._k):
                    break

            self.print_verbose('')
            self._k += 1
        return True

    def run(self, *args, **kwargs):
        self.run_initialize(*args, **kwargs)
        return self.run_loop()
