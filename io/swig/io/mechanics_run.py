

# Mechanics IO

"""Run a pre-configured Siconos "mechanics" HDF5 file."""

from __future__ import print_function

import os
import sys

from math import cos, sin, asin, atan2
from scipy import constants

import numpy as np
import h5py
import bisect
import time
import numbers
import shutil


import tempfile
from contextlib import contextmanager

# Siconos imports
import siconos.io.mechanics_hdf5
import siconos.numerics as Numerics
from siconos.kernel import \
    EqualityConditionNSL, \
    Interaction, DynamicalSystem, TimeStepping,\
    SICONOS_OSNSP_TS_VELOCITY,\
    cast_FrictionContact
import siconos.kernel as Kernel

# Siconos Mechanics imports
from siconos.mechanics.collision.tools import Contactor, Shape
from siconos.mechanics import joints
from siconos.io.io_base import MechanicsIO

# Imports for mechanics 'collision' submodule
from siconos.mechanics.collision import RigidBodyDS, RigidBody2dDS, \
    SiconosSphere, SiconosBox, SiconosCylinder, SiconosCone, SiconosCapsule, SiconosPlane,  SiconosConvexHull, \
    SiconosDisk, SiconosBox2d, SiconosConvexHull2d, \
    SiconosContactor, SiconosContactorSet, \
    SiconosMesh, SiconosHeightMap

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
    default_simulation_class = TimeStepping
    default_body_class = RigidBodyDS
    if backend == 'bullet':
        if not have_bullet:
            print('[mechanics_run] WARNING : The backend for collision is bullet but it seems that bullet is not linked with siconos')
        def m(options):
            if options is None:
                options = SiconosBulletOptions()
            return SiconosBulletCollisionManager(options)
        default_manager_class = m
        use_bullet = have_bullet
    elif backend == 'occ':
        default_manager_class = lambda options: occ.OccSpaceFilter()
        default_simulation_class = occ.OccTimeStepping
        default_body_class = occ.OccBody
        use_bullet = False

setup_default_classes()

### Constants

joint_points_axes = {
    'KneeJointR': (1, 0),
    'PivotJointR': (1, 1),
    'PrismaticJointR': (0, 1),
    'CylindricalJointR': (1, 1),
    'FixedJointR': (0, 0),
}

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


class Timer():

    def __init__(self):
        self._t0 = time.clock()

    def elapsed(self):
        return time.clock() - self._t0

    def update(self):
        self._t0 = time.clock()


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
                                chunks=[None,(4000,nbcolumns)][comp],
                                compression=[None,'gzip'][comp],
                                compression_opts=[None,9][comp])


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
            if exc.errno != errno.EEXIST:
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
            axis=orientation[0]
            if not (len(axis) == 3):
                raise AssertionError("quaternion_get. axis in not a 3D vector")
            angle=orientation[1]
            if not (isinstance(angle,float)):
                raise AssertionError("quaternion_get. angle must be a float")
            n=sin(angle / 2.) / np.linalg.norm(axis)
            ori=[cos(angle / 2.), axis[0] * n, axis[1] * n, axis[2] * n]
    else:
        if not (len(orientation)==4):
            raise AssertionError("quaternion_get. the quaternion must be of size 4")
        # a given quaternion
        ori=orientation
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
        self._primitive = {'Sphere': SiconosSphere,
                           'Box': SiconosBox,
                           'Cylinder': SiconosCylinder,
                           'Capsule': SiconosCapsule,
                           'Cone':SiconosCone,
                           'Plane': SiconosPlane,
                           'Disk': SiconosDisk,
                           'Box2d':SiconosBox2d}

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
                            scale = self.attributes(shape_name).get('scale',None)
                            mesh, dims = load_siconos_mesh(tmpf[1], scale=scale)
                            self._shapes[shape_name] = mesh
                            mesh.setInsideMargin(
                                self.shape(shape_name).attrs.get('insideMargin',
                                                                 min(dims)*0.02))
                            mesh.setOutsideMargin(
                                self.shape(shape_name).attrs.get('outsideMargin',0))
                    else:
                        raise AssertionError("ShapeCollection.get() self.attributes(shape_name)['type'] != 'vtp'")

                elif self.attributes(shape_name)['type'] in ['step', 'stp']:
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



                            #ok = step_reader.TransferRoot(1)
                            step_reader.TransferRoots() # VA : We decide to loads all shapes in the step file
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

                    if not (self.shape(shape_name).dtype == h5py.new_vlen(str)):
                        raise AssertionError("ShapeCollection.get()")

                    #assert(self.shape(shape_name).dtype == h5py.new_vlen(str))

                    with tmpfile(contents=self.shape(shape_name)[:][0]) as tmpf:
                        iges_reader = IGESControl_Reader()

                        status = iges_reader.ReadFile(tmpf[1])

                        if status == IFSelect_RetDone:  # check status
                            failsonly = False
                            iges_reader.PrintCheckLoad(
                                failsonly, IFSelect_ItemsByEntity)
                            iges_reader.PrintCheckTransfer(
                                failsonly, IFSelect_ItemsByEntity)

                            ok = iges_reader.TransferRoots()
                            nbs = iges_reader.NbShapes()

                            for i in range(1, nbs + 1):
                                shape = iges_reader.Shape(i)
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
                    # a heightmap defined by a matrix
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

                elif self.attributes(shape_name)['type'] in ['convex']:
                    # a convex point set
                    points = self.shape(shape_name)
                    if self._io._dimension == 3:
                        convex = SiconosConvexHull(points)
                        dims = [points[:,0].max() - points[:,0].min(),
                                points[:,1].max() - points[:,1].min(),
                                points[:,2].max() - points[:,2].min()]
                    elif self._io._dimension == 2:
                        convex = SiconosConvexHull2d(points)
                        dims = [points[:,0].max() - points[:,0].min(),
                                points[:,1].max() - points[:,1].min()]
                        
                    convex.setInsideMargin(
                        self.shape(shape_name).attrs.get('insideMargin',
                                                         min(dims)*0.02))
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


class MechanicsHdf5Runner(siconos.io.mechanics_hdf5.MechanicsHdf5):

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
                 interaction_manager=None, nsds=None, simulation=None,
                 osi=None, shape_filename=None,
                 set_external_forces=None, gravity_scale=None, collision_margin=None,
                 use_compression=False, output_domains=False, verbose=True):
        
        super(MechanicsHdf5Runner, self).__init__(io_filename, mode, None,
                                                  use_compression, output_domains, verbose)
        self._interman = interaction_manager
        self._nsds = nsds
        self._simulation = simulation
        self._osi = osi
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
        self._contact_index_set = 1
        self._scheduled_births = []
        self._start_run_iteration_hook = None
        self._end_run_iteration_hook = None

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
        if self._dimension == 3:
            weight = [0, 0, - body.scalarMass() * g]
        elif self._dimension == 2:
            weight = [0, - body.scalarMass() * g, 0.]
        body.setFExtPtr(weight)

    def import_nonsmooth_law(self, name):
        if self._interman is not None:
            nslawClass = getattr(Kernel, self._nslaws_data[name].attrs['type'])
            if nslawClass == Kernel.NewtonImpactFrictionNSL:
                nslaw = nslawClass(float(self._nslaws_data[name].attrs['e']), 0.,
                                   float(self._nslaws_data[name].attrs['mu']), self._dimension)
            elif nslawClass == Kernel.NewtonImpactRollingFrictionNSL:
                nslaw = nslawClass(float(self._nslaws_data[name].attrs['e']), 0.,
                                   float(self._nslaws_data[name].attrs['mu']),
                                   float(self._nslaws_data[name].attrs['mu_r']), 5)
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
                self._interman.insertNonSmoothLaw(nslaw, gid1, gid2)

    def import_occ_object(self, name, translation, orientation,
                        velocity, contactors, mass, given_inertia, body_class,
                        shape_class, face_class, edge_class, birth=False, number=None):

        if mass is None :
            # a static object
            body = None

            self._static[name] = {
                    'number': number,
                    'origin': translation,
                    'orientation': orientation}
        else:
            if body_class is None:
                body_class = occ.OccBody

            assert (given_inertia is not None)
            inertia = given_inertia
            if inertia is not None:
                if np.shape(inertia) == (3,):
                    inertia = np.diag(inertia)
                elif np.shape(inertia) != (3,3):
                    self.print_verbose('Wrong shape of inertia')


            body = body_class(
                list(translation) + list(orientation), velocity, mass, inertia)

            if number is not None:
                body.setNumber(number)

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
            self._nsds.insertDynamicalSystem(body)
            self._nsds.setName(body, str(name))

            if birth and self._verbose:
                self.print_verbose ('birth of body named {0}, translation {1}, orientation {2}'.format(name, translation, orientation))

        return body

    def import_bullet_object(self, name, translation, orientation,
                           velocity, contactors, mass, inertia,
                           body_class, shape_class, birth=False,
                           number = None):


        if body_class is None:
            body_class = default_body_class

        if self._dimension ==2:
            body_class = RigidBody2dDS
            
        if self._interman is not None and 'input' in self._data:
            body = None
            if mass is None :
                # a static object
                cset = SiconosContactorSet()
                csetpos = (translation + orientation)
                for c in contactors:
                    shp = self._shape.get(c.shape_name)
                    pos = list(c.translation) + list(c.orientation)
                    cset.append(SiconosContactor(shp, pos, c.group))
                    self.print_verbose('Adding shape %s to static contactor'%c.shape_name, pos)

                self._interman.insertStaticContactorSet(cset, csetpos)

                self._static[name] = {
                    'number': number,
                    'origin': translation,
                    'orientation': orientation,
                    'shape': shp,
                }

            else: # dynamic object

                if not np.isscalar(mass) or  mass <=0:
                    self.print_verbose('**** Warning mass must be positive scalar')

                inertia_ok=False
                if self._dimension ==3:
                    if inertia is not None:
                        if np.shape(inertia) == (3,):
                            inertia = np.diag(inertia)
                            inertia_ok = True
                        if np.shape(inertia) == (3,3):
                            inertia_ok = True
                elif self._dimension ==2:
                    if inertia is not None:
                        if np.shape(inertia) == (1,1) or  np.shape(inertia) == (1,) or np.isscalar(inertia) :
                            inertia_ok = True


                if inertia_ok:
                    # create the dynamics object with mass and inertia
                    body = body_class(translation + orientation,
                                      velocity, mass, inertia)
                    body.setUseContactorInertia(False)
                else:
                    if inertia is not None:
                        if self._dimension ==3:
                            self.print_verbose('**** Warning inertia for object named {0} does not have the correct shape: {1} instead of (3, 3) or (3,)'.format(name, np.shape(inertia)))
                            self.print_verbose('**** Inertia will be computed with the shape of the first contactor')
                        elif self._dimension ==2:
                            self.print_verbose('**** Warning inertia for object named {0} does not have the correct shape: {1} instead of (1, 1) or (1,) or scalar'.format(name, np.shape(inertia)))
                            self.print_verbose('**** Inertia will be computed with the shape of the first contactor')
                            
                    body = body_class(translation + orientation,
                                      velocity, mass)
                    body.setUseContactorInertia(True)


                self_collide = self._input[name].get('allow_self_collide',None)
                if self_collide is not None:
                    body.setAllowSelfCollide(not not self_collide)

                cset = SiconosContactorSet()
                for c in contactors:
                    shp = self._shape.get(c.shape_name)
                    pos = list(c.translation) + list(c.orientation)
                    cset.append(SiconosContactor(shp, pos, c.group))

                body.setContactors(cset)

            if body:
                # set id number
                if number is not None:
                    body.setNumber(number)

                # set external forces
                self._set_external_forces(body)

                # add the dynamical system to the non smooth
                # dynamical system
                self._nsds.insertDynamicalSystem(body)
                self._nsds.setName(body, str(name))

        return body

    def make_coupler_jointr(self, ds1_name, ds2_name, coupled, references):
        topo = self._nsds.topology()
        dof1, dof2, ratio = coupled[0,:]
        refds_name = None
        if len(references)>0:
            joint1_name = references[0]
            joint1 = joints.cast_NewtonEulerJointR(
                topo.getInteraction(str(joint1_name)).relation())
            joint2 = joint1
            joint1_ds1 = self.joints()[joint1_name].attrs['object1']
            joint1_ds2 = self.joints()[joint1_name].attrs.get('object2',None)
        if len(references)>1:
            joint2_name = references[1]
            joint2 = joints.cast_NewtonEulerJointR(
                topo.getInteraction(str(joint2_name)).relation())
            joint2_ds1 = self.joints()[joint2_name].attrs['object1']
            joint2_ds2 = self.joints()[joint2_name].attrs.get('object2',None)

            if len(references)>2:
                refds_name = references[2]
            else:
                # if second joint provided but no reference, then
                # infer refds_name from the reference joints
                dss = set([joint1_ds1, joint1_ds2, joint2_ds1, joint2_ds2])
                diff = list(dss.difference(set([ds1_name, ds2_name])))
                # there must be exactly one reference in common that
                # is not either of the DSs
                assert(len(diff)==1)
                refds_name = diff[0]

        if refds_name:
            refds = topo.getDynamicalSystem(str(refds_name))

            # Determine reference indexes:
            # Assert if neither ds in reference joints is the
            # ref ds and that the other ds is ds1/ds2 as
            # appropriate.
            refidx1, refidx2 = 0, 0
            if joint1_ds1 == ds1_name:
                refidx1 = 1
                assert(joint1_ds2 == refds_name)
            else:
                assert(joint1_ds1 == refds_name)
                assert(joint1_ds2 == ds1_name)
            if joint2_ds1 == ds2_name:
                refidx2 = 1
                assert(joint2_ds2 == refds_name)
            else:
                assert(joint2_ds1 == refds_name)
                assert(joint2_ds2 == ds2_name)

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
                'allow_self_collide',None)
            stops = self.joints()[name].attrs.get('stops',None)
            nslaws = self.joints()[name].attrs.get('nslaws',None)
            friction = self.joints()[name].attrs.get('friction',None)
            coupled = self.joints()[name].attrs.get('coupled',None)
            references = self.joints()[name].attrs.get('references',None)

            # work-around h5py unicode bug
            # https://github.com/h5py/h5py/issues/379
            if references is not None:
                references = [r.decode('utf-8') if isinstance(r,bytes) else r
                              for r in references]

            points = self.joints()[name].attrs.get('points',[])
            axes = self.joints()[name].attrs.get('axes',[])
            siconos.io.mechanics_hdf5.check_points_axes(name, joint_class, points, axes)

            ds1_name = self.joints()[name].attrs['object1']
            ds1 = topo.getDynamicalSystem(ds1_name)
            ds2 = None

            if 'object2' in self.joints()[name].attrs:
                ds2_name = self.joints()[name].attrs['object2']
                ds2 = topo.getDynamicalSystem(ds2_name)

            if joint_class == joints.CouplerJointR:
                # This case is a little different, handle it specially
                assert(references is not None)
                assert(np.shape(coupled)==(1,3))
                joint = self.make_coupler_jointr(ds1_name, ds2_name,
                                                coupled, references)
                coupled = None # Skip the use for "coupled" below, to
                               # install joint-local couplings

            else:
                # Generic NewtonEulerJointR interface
                joint = joint_class()
                for n,p in enumerate(points):
                    joint.setPoint(n, p)
                for n,a in enumerate(axes):
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
                assert np.shape(stops)[1] == 3, 'Joint stops shape must be (?,3)'
                if nslaws is None:
                    nslaws = [Kernel.NewtonImpactNSL(0.0)]*np.shape(stops)[0]
                elif isinstance(nslaws,bytes):
                    nslaws = [self._nslaws[nslaws.decode('utf-8')]]*np.shape(stops)[0]
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
                    self._nsds.\
                        link(stop_inter, ds1, ds2)
                    nsds.setName(stop_inter, '%s_stop%d'%(str(name),n))

            # The per-axis friction NSL, can be ''
            if friction is not None:
                if isinstance(friction,str):
                    friction = [friction]
                elif isinstance(friction,bytes):
                    friction = [friction.decode('utf-8')]
                else:
                    friction = [(f.decode('utf-8') if isinstance(f,bytes) else f)
                                for f in friction]
                for ax,fr_nslaw in enumerate(friction):
                    if fr_nslaw == '':  # no None in hdf5, use empty string
                        continue        # instead for no NSL on an axis
                    nslaw = self._nslaws[fr_nslaw]
                    fr = joints.JointFrictionR(joint, [ax])
                    fr_inter = Interaction(nslaw, fr)
                    self._nsds.\
                        link(fr_inter, ds1, ds2)
                    nsds.setName(fr_inter, '%s_friction%d'%(str(name),ax))

            # An array of tuples (dof1, dof2, ratio) specifies
            # coupling between a joint's DoFs (e.g., to turn a
            # cylindrical joint into a screw joint)
            if coupled is not None:
                if len(coupled.shape)==1: coupled = numpy.array([coupled])
                assert(coupled.shape[1]==3)
                for n, (dof1, dof2, ratio) in enumerate(coupled):
                    cpl = joints.CouplerJointR(joint, int(dof1),
                                               joint, int(dof2), ratio)
                    cpl.setBasePositions(q1, q2)
                    cpl_inter = Interaction(EqualityConditionNSL(1), cpl)
                    self._nsds.\
                        link(cpl_inter, ds1, ds2)
                    nsds.setName(cpl_inter, '%s_coupler%d'%(str(name),n))

    def import_boundary_conditions(self, name):
        if self._interman is not None:
            topo = self._nsds.\
                topology()

            bc_type = self.boundary_conditions()[name].attrs['type']
            bc_class = getattr(Kernel,bc_type)

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
            #    self._nsds.\
            #        link(joint_inter, ds1)

    def import_permanent_interactions(self, name):
        """
        """
        if (self._interman is not None and 'input' in self._data
              and self.permanent_interactions() is not None):
            topo = self._nsds.\
                topology()

            pinter = self.permanent_interactions()[name]
            body1_name = pinter.attrs['body1_name']
            body2_name = pinter.attrs['body2_name']

            try:
                ds1 = \
                      Kernel.cast_NewtonEulerDS(
                          topo.getDynamicalSystem(body1_name))
            except:
                ds1 = None

            try:
                ds2 = \
                      Kernel.cast_NewtonEulerDS(
                          topo.getDynamicalSystem(body2_name))
            except Exception :
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
                        self._nsds.link(inter, ds1,
                                                                    ds2)
                    else:
                        self._nsds.link(inter, ds1)

                    # keep pointers
                    self._keep.append([cocs1, cocs2, cp1,
                                       cp2, relation])

    def import_object(self, name, body_class=None, shape_class=None,
                     face_class=None, edge_class=None, birth=False,
                     translation=None, orientation=None, velocity=None):
        """
        Import an object by name, possibly overriding initial position and velocity.
        """
        obj = self._input[name]
        self.print_verbose ('Import  dynamic or static object number ',
                            obj.attrs['id'], 'from initial state')
        self.print_verbose ('                object name   ', name)


        if translation is None:
            translation = obj.attrs['translation']
        if orientation is None:
            orientation = obj.attrs['orientation']
        if velocity is None:
            velocity = obj.attrs['velocity']

        # bodyframe center of mass
        center_of_mass = floatv(obj.attrs.get('center_of_mass', [0,0,0]))

        mass = obj.attrs.get('mass', None)
        inertia = obj.attrs.get('inertia', None)

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
                        relative_translation=np.subtract(ctr.attrs['translation'].astype(float), center_of_mass),
                        relative_orientation=ctr.attrs['orientation'].astype(float)))
            elif 'group' in ctr.attrs:
                # bullet contact
                assert not occ_type
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
            body = self.import_occ_object(
                name, floatv(translation), floatv(orientation),
                floatv(velocity), contactors, mass,
                inertia, body_class, shape_class, face_class,
                edge_class, birth=birth,
                number = self.instances()[name].attrs['id'])
        else:
            # Bullet object
            body = self.import_bullet_object(
                name, floatv(translation), floatv(orientation),
                floatv(velocity), contactors, mass,
                inertia, body_class, shape_class, birth=birth,
                number = self.instances()[name].attrs['id'])

        # schedule its death immediately
        time_of_death = obj.attrs.get('time_of_death', None)
        if None not in (time_of_death, body):
            bisect.insort_left(self._scheduled_deaths, time_of_death)
            if time_of_death in self._deaths:
                self._deaths[time_of_death].append((name, obj, body))
            else:
                self._deaths[time_of_death] = [(name, obj, body)]

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



        
        for shape_name in self._ref:
            self._number_of_shapes += 1

        # import dynamical systems
        if self._interman is not None and 'input' in self._data:

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
                id_vlast = np.where(
                            abs(velocities[:, 0] - max_time) < 1e-9)[0]

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

            for (name, obj) in sorted(self._input.items(),
                                      key=lambda x: x[0]):

                mass = obj.attrs.get('mass', None)
                time_of_birth = obj.attrs.get('time_of_birth',-1)

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
                    if (mass is not None
                        and dpos_data is not None and len(dpos_data) > 0):

                        self.print_verbose('Import  dynamic object name ',
                                           name,
                                           'from current state')
                        self.print_verbose('imported object has id: {0}'.
                                           format(obj.attrs['id']))

                        xpos = xdpos_data[obj.attrs['id']]
                        translation = (xpos[2], xpos[3], xpos[4])
                        orientation = (xpos[5], xpos[6], xpos[7], xpos[8])

                        xvel = xvelocities[obj.attrs['id']]
                        velocity = (xvel[2], xvel[3], xvel[4],
                                    xvel[5], xvel[6], xvel[7])

                        self.print_verbose('position:', list(translation) +
                                           list(orientation))
                        self.print_verbose('velocity:',  velocity)


                    else:
                        # start from initial conditions
                        self.print_verbose ('Import  dynamic or static object number ', obj.attrs['id'], 'from initial state')
                        self.print_verbose ('                object name   ', name)
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
            for name in self._nslaws_data:
                self.import_nonsmooth_law(name)

            for name in self.joints():
                self.import_joint(name)

            for name in self.boundary_conditions():
                self.import_boundary_conditions(name)

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
            for (name, obj) in self._births[time_of_birth]:
                self.import_object(name, body_class, shape_class,
                                  face_class, edge_class, birth=True)

    def execute_deaths(self):
        """
        Remove objects from the NSDS
        """
        time = self.current_time()

        ind_time = bisect.bisect_right(self._scheduled_deaths, time)

        current_times_of_deaths = set(self._scheduled_deaths[:ind_time])
        self._scheduled_deaths = self._scheduled_deaths[ind_time:]

        for time_of_death in current_times_of_deaths:
            for (name, obj, body) in self._deaths[time_of_death]:
                self._interman.removeBody(body)
                self._nsds.removeDynamicalSystem(body)

    def output_static_objects(self):
        """
        Outputs translations and orientations of static objects
        """
        time = self.current_time()
        p = 0
        self._static_data.resize(len(self._static), 0)

        for static in self._static.values():

            self.print_verbose('output static object', static['number'])

            self.print_verbose (static.keys())
            translation = static['origin']
            rotation = static['orientation']


            if self._dimension ==3 :
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
           
            elif self._dimension ==2 :
                # VA. change the position such that is corresponds to a 3D object
                self._static_data[p, :] = \
                    [time,
                     static['number'],
                     translation[0],
                     translation[1],
                     0.0,
                     cos(rotation[0]/2.0),
                     0.0,
                     0.0,
                     sin(rotation[0]/2.0)]
            
            p += 1

    def output_dynamic_objects(self, initial=False):
        """
        Outputs translations and orientations of dynamic objects.
        """

        current_line = self._dynamic_data.shape[0]

        time = self.current_time()

        positions = self._io.positions(self._nsds)

        if positions is not None:
            
            self._dynamic_data.resize(current_line + positions.shape[0], 0)

            times = np.empty((positions.shape[0], 1))
            times.fill(time)
            
            if self._dimension == 3 :
                self._dynamic_data[current_line:, :] = np.concatenate((times,
                                                                       positions),
                                                                      axis=1)
            elif self._dimension == 2 :
                 # VA. change the position such that is corresponds to a 3D object
                new_positions=np.zeros((positions.shape[0],8))

                new_positions[:,0]=positions[:,0]
                new_positions[:,1]=positions[:,1]
                new_positions[:,2]=positions[:,2]

                new_positions[:,4]=np.cos(positions[:,3]/2.0)
                new_positions[:,7]=np.sin(positions[:,3]/2.0)
                self._dynamic_data[current_line:, :] = np.concatenate((times,
                                                                       new_positions),
                                                                      axis=1)

                

    def output_velocities(self):
        """
        Output velocities of dynamic objects
        """

        current_line = self._dynamic_data.shape[0]

        time = self.current_time()

        velocities = self._io.velocities(self._nsds)

        if velocities is not None:

            self._velocities_data.resize(current_line + velocities.shape[0], 0)

            times = np.empty((velocities.shape[0], 1))
            times.fill(time)
            if self._dimension == 3 :
                self._velocities_data[current_line:, :] = np.concatenate((times,
                                                                          velocities),
                                                                         axis=1)
            elif self._dimension == 2 :
                 # VA. change the position such that is corresponds to a 3D object
                new_velocities=np.zeros((velocities.shape[0],7))

                new_velocities[:,0]=velocities[:,0]
                
                new_velocities[:,1]=velocities[:,1]
                new_velocities[:,2]=velocities[:,2]

                new_velocities[:,6]=velocities[:,3]
                self._dynamic_data[current_line:, :] = np.concatenate((times,
                                                                       new_velocities),
                                                                      axis=1)

                
    def output_contact_forces(self):
        """
        Outputs contact forces
        _contact_index_set default value is 1.
        """
        if self._nsds.\
                topology().indexSetsSize() > 1:
            time = self.current_time()
            contact_points = self._io.contactPoints(self._nsds,
                                                    self._contact_index_set)

            if contact_points is not None:
                current_line = self._cf_data.shape[0]
                # Increase the number of lines in cf_data
                # (h5 dataset with chunks)
                self._cf_data.resize(current_line + contact_points.shape[0], 0)
                times = np.empty((contact_points.shape[0], 1))
                times.fill(time)

                if self._dimension == 3:
                    new_contact_points = np.zeros((contact_points.shape[0],25))
                     
                    new_contact_points[:, 0] = contact_points[:,0] # mu
                    new_contact_points[:, 1] = contact_points[:,1] # posa
                    new_contact_points[:, 2] = contact_points[:,2]
                    #new_contact_points[:,3]
                    new_contact_points[:,4] = contact_points[:,3] # posb
                    new_contact_points[:,5] = contact_points[:,4]
                    #new_contact_points[:,6]
                    new_contact_points[:,7] = contact_points[:,5] # nc
                    new_contact_points[:,8] = contact_points[:,6]
                    #new_contact_points[:,9]
                    new_contact_points[:,10] = contact_points[:,7] # cf 
                    new_contact_points[:,11] = contact_points[:,8]
                    #new_contact_points[:,12]
                    new_contact_points[:,13] = contact_points[:,9] # gap
                    new_contact_points[:,14] = contact_points[:,10]
                    #new_contact_points[:,15]
                    new_contact_points[:,16] = contact_points[:,11] # relative velocity
                    new_contact_points[:,17] = contact_points[:,12]
                    #new_contact_points[:,18]
                    new_contact_points[:,19] = contact_points[:,13] # reaction impulse
                    new_contact_points[:,20] = contact_points[:,14]
                    #new_contact_points[:,21]
                    new_contact_points[:,22] = contact_points[:,15] # inter id
                    new_contact_points[:,23] = contact_points[:,16] # ds 1
                    new_contact_points[:,24] = contact_points[:,17] # ds 2
        
                    self._cf_data[current_line:, :] = \
                        np.concatenate((times,
                                        contact_points),
                                       axis=1)

                elif self._dimension == 2:
                    
                    # VA. change the contact info such that is corresponds to a 3D object
                     new_contact_points=np.zeros((contact_points.shape[0],25))
                     
                     new_contact_points[:,0] = contact_points[:,0] # mu
                     new_contact_points[:,1] = contact_points[:,1] # posa
                     new_contact_points[:,2] = contact_points[:,2]
                     #new_contact_points[:,3]
                     new_contact_points[:,4] = contact_points[:,3] # posb
                     new_contact_points[:,5] = contact_points[:,4]
                     #new_contact_points[:,6]
                     new_contact_points[:,7] = contact_points[:,5] # nc
                     new_contact_points[:,8] = contact_points[:,6]
                     #new_contact_points[:,9]
                     new_contact_points[:,10] = contact_points[:,7] # cf 
                     new_contact_points[:,11] = contact_points[:,8]
                     #new_contact_points[:,12]
                     new_contact_points[:,13] = contact_points[:,9] # gap
                     new_contact_points[:,14] = contact_points[:,10]
                     #new_contact_points[:,15]
                     new_contact_points[:,16] = contact_points[:,11] # relative velocity
                     new_contact_points[:,17] = contact_points[:,12]
                     #new_contact_points[:,18]
                     new_contact_points[:,19] = contact_points[:,13] # reaction impulse
                     new_contact_points[:,20] = contact_points[:,14]
                     #new_contact_points[:,21]
                     new_contact_points[:,22] = contact_points[:,15] # inter id
                     new_contact_points[:,23] = contact_points[:,16] # ds 1
                     new_contact_points[:,24] = contact_points[:,17] # ds 2
                     self._cf_data[current_line:, :] = \
                         np.concatenate((times,
                                         new_contact_points),
                                        axis=1)
                     

                # return the number of contacts
                return len(contact_points)

            


            
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
        so = self._simulation.oneStepNSProblem(0).\
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

    def print_solver_infos(self):
        """
        Outputs solver #iterations & precision reached
        """
        time = self.current_time()
        so = self._simulation.oneStepNSProblem(0).\
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


        self.print_verbose('Numerics solver info at time : {0:10.6f}'.format(time),
              'iterations = {0:8d}'.format(iterations),
              'precision = {0:5.3e}'.format(precision))

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
        subprocess.check_call(['siconos', '--noexec','.'])

    def import_external_functions(self):
        topo = self._nsds.\
                topology()

        for name in self._external_functions:

            ext_fun = self._external_functions[name]
            plugin_name = ext_fun.attrs['plugin_name']
            plugin_function_name = ext_fun.attrs['plugin_function_name']
            body_name = ext_fun.attrs['body_name']
            ds = Kernel.cast_NewtonEulerDS(topo.getDynamicalSystem(body_name))

            if 'function_name' in ext_fun.attrs:
                function_name = ext_fun.attrs['function_name']

                getattr(ds, function_name)(plugin_name,
                                           plugin_function_name)
            else:
                bc_indices = ext_fun.attrs['bc_indices']
                # a boundary condition
                bc = Kernel.BoundaryCondition(bc_indices)
                bc.setComputePrescribedVelocityFunction(plugin_name,
                                                        plugin_function_name)
                ds.setBoundaryConditions(bc)

    def explode_Newton_solve(self,  with_timer):
        s = self._simulation
        self.log(s.initialize, with_timer)()
        self.log(s.resetLambdas, with_timer)()
        # Again the access to s._newtonTolerance generates a segfault due to director
        newtonTolerance = s.newtonTolerance()
        newtonMaxIteration = s.newtonMaxIteration()

        # return _kernel.TimeStepping_newtonSolve(self, criterion, maxStep)
        # RuntimeError: accessing protected member newtonSolve
        #s.newtonSolve(newtonTolerance, newtonMaxIteration);

        newtonNbIterations =0
        isNewtonConverge =False
        explode_computeOneStep = False

        self.log(s.initializeNewtonLoop, with_timer)()
        while (not isNewtonConverge) and (newtonNbIterations <newtonMaxIteration):
            #self.print_verbose('newtonNbIterations',newtonNbIterations)
            info=0
            newtonNbIterations = newtonNbIterations+1
            self.log(s.prepareNewtonIteration, with_timer)()
            self.log(s.computeFreeState, with_timer)()
            if s.numberOfOSNSProblems() >0:
                if explode_computeOneStep:
                    # experimental
                    osnsp=s.oneStepNSProblem(SICONOS_OSNSP_TS_VELOCITY)
                    #info = self.log(osnsp.compute, with_timer)(s.nextTime())
                    fc = cast_FrictionContact(osnsp)
                    #self.log(fc.updateInteractionBlocks, with_timer)()
                    self.log(fc.preCompute, with_timer)(s.nextTime())
                    self.log(fc.updateMu, with_timer)()
                    info = self.log(fc.solve, with_timer)()
                    self.log(fc.postCompute, with_timer)()
                else:
                    info = self.log(s.computeOneStepNSProblem, with_timer)(SICONOS_OSNSP_TS_VELOCITY)
            self.log(s.DefaultCheckSolverOutput, with_timer)(info);
            self.log(s.updateInput, with_timer)();
            self.log(s.updateState, with_timer)();
            if (not isNewtonConverge) and (newtonNbIterations <newtonMaxIteration):
                self.log(s.updateOutput, with_timer)()
            isNewtonConverge = self.log(s.newtonCheckConvergence, with_timer)
            (newtonTolerance)
            if s.displayNewtonConvergence():
                s.displayNewtonConvergenceInTheLoop();
            if (not isNewtonConverge) and (not info):
                if s.numberOfOSNSProblems() > 0:
                    self.log(s.saveYandLambdaInOldVariables, with_timer)()
        if s.displayNewtonConvergence():
            s.displayNewtonConvergenceAtTheEnd(info, newtonMaxIteration);



    def run(self,
            with_timer=False,
            time_stepping=None,
            interaction_manager=None,
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
            local_solver=Numerics.SICONOS_FRICTION_3D_ONECONTACT_NSN_GP_HYBRID,
            itermax=100000,
            osnspb_max_size=0,
            tolerance=1e-8,
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
            friction_contact_trace=False,
            friction_contact_trace_params=None,
            contact_index_set=1,
            osi=Kernel.MoreauJeanOSI,
            constraint_activation_threshold=0.0,
            explode_Newton_solve=False,
            display_Newton_convergence=False,
            start_run_iteration_hook=None,
            end_run_iteration_hook=None
            ):
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
          multiPoint_iterations : use bullet "multipoint iterations"
                                 (default True)
          theta : parameter for Moreau-Jean OSI (default 0.50001)
          Newton_max_iter : maximum number of iterations for
                          integrator Newton loop (default 20)
          set_external_forces : method for external forces
                                (default earth gravity)
          solver : OneStepNsProblem solver  (default Numerics.SICONOS_FRICTION_3D_NSGS)
          local_solver : OneStepNsProblem solver local solver (default Numerics.SICONOS_FRICTION_3D_ONECONTACT_NSN_GP_HYBRID)
          itermax : maximum number of iteration for solver (default 100000)
          osnspb_max_size : estimation of maximum number of dynamical systems for memory pre-allocations.
                            if 0, set to simulation().nonSmoothDynamicalSystem().topology().numberOfConstraints()
                            (default 0)
          tolerance : friction contact solver tolerance (default 1e-8)
          exit_tolerance : if not None, the simulation will stop if precision >= exit_tolerance (default None)
          numerics_verbose : set verbose mode in numerics
          output_frequency : 0 to disable (default 1)
          contact_index_set : index set from which contact point information is retrieved.
        """
        self._verbose = verbose
        self.print_verbose ('load siconos module ...')
        from siconos.kernel import \
            NonSmoothDynamicalSystem,\
            TimeDiscretisation,\
            GenericMechanical, FrictionContact,\
            GlobalFrictionContact, RollingFrictionContact,\
            NewtonImpactFrictionNSL, NewtonImpactRollingFrictionNSL

        from siconos.numerics import SICONOS_FRICTION_3D_ONECONTACT_NSN
        from siconos.numerics import SICONOS_GLOBAL_FRICTION_3D_ADMM

        self.print_verbose ('setup model simulation ...')
        if set_external_forces is not None:
            self._set_external_forces=set_external_forces

        if interaction_manager is None:  interaction_manager = default_manager_class
        if time_stepping is None: time_stepping = default_simulation_class

        if output_frequency is not None:
            self._output_frequency=output_frequency

        if output_backup_frequency is not None:
            self._output_backup_frequency=output_backup_frequency

        if output_backup is not None:
            self._output_backup=output_backup
 
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
            self.print_verbose('')
            self.print_verbose('Simulation will run from {0:.4f} to {1:.4f}s, step {2} to step {3} (h={4}, times=[{5},{6}])'
                                       .format(t0, T, k0, kT, h,
                                               min(times) if len(times)>0 else '?',
                                               max(times) if len(times)>0 else '?'))
            self.print_verbose('')
        else:
            self.print_verbose('Simulation time {0} >= T={1}, exiting.'.format(t0,T))
            return

        # get dimension
        self._dimension=self._out.attrs.get('dimension', 3)

        # Respect run() parameter for multipoints_iterations for
        # backwards compatibility, but this is overridden by
        # SiconosBulletOptions if one is provided.
        if multipoints_iterations is not None and options is None:
            options = SiconosBulletOptions()
            options.perturbationIterations = 3*multipoints_iterations
            options.minimumPointsPerturbationThreshold = 3*multipoints_iterations

        if (self._dimension ==2) :
            if  options is None:
                options = SiconosBulletOptions()
            options.dimension = SICONOS_BULLET_2D
            
        self._interman = interaction_manager(options)

        joints = list(self.joints())
        if hasattr(self._interman, 'useEqualityConstraints') and len(joints)==0:
            self._interman.useEqualityConstraints(False)

        # (0) NonSmooth Dynamical Systems definition
        self._nsds=NonSmoothDynamicalSystem(t0, T)
        nsds=self._nsds

        self.print_verbose ('import scene ...')
        self.import_scene(t0, body_class, shape_class, face_class, edge_class)

        self._contact_index_set = contact_index_set

        # (1) OneStepIntegrators


        self._osi=osi(theta)
        if (osi == Kernel.MoreauJeanOSI):
            self._osi.setConstraintActivationThreshold(constraint_activation_threshold)
        # (2) Time discretisation --
        timedisc=TimeDiscretisation(t0, h)


        # (3) choice of default OneStepNonSmoothProblem w.r.t the type of nslaws
        nslaw_type_list =[]
        for name in self._nslaws_data:
            nslaw_type_list.append(self._nslaws_data[name].attrs['type'])

        #print(set(nslaw_type_list))

        # This trick is used to add the EqualityConditionNSL to the list of nslaw type
        # this must be improved by adding the EqualityConditionNSL in self._nslaws_data
        # when a joint is imported.
        # For the moment, the nslaw is implicitely added when we import_joint but is not stored
        # self._nslaws_data

        if len(joints) > 0:
            nslaw_type_list.append('EqualityConditionNSL')


        nb_of_nslaw_type =  len(set(nslaw_type_list))
        # print(set(nslaw_type_list))
        # input()
        if (friction_contact_trace == False) :
            if (osi == Kernel.MoreauJeanGOSI):
                if (nb_of_nslaw_type >1) or 'NewtonImpactFrictionNSL' not in set(nslaw_type_list):
                    raise RuntimeError("MoreauJeanGOSI can deal only with NewtonImpactFrictionNSL nslaw")
                else:
                    if (solver == Numerics.SICONOS_FRICTION_3D_NSGS):
                        osnspb=GlobalFrictionContact(3,SICONOS_GLOBAL_FRICTION_3D_ADMM)
                        osnspb.setMStorageType(2)
                    else:
                        osnspb=GlobalFrictionContact(3,solver)
                        osnspb.setMStorageType(1)
                    osnspb.setMaxSize(osnspb_max_size)
            else:
                if (nb_of_nslaw_type >1):
                    osnspb=GenericMechanical(SICONOS_FRICTION_3D_ONECONTACT_NSN)
                    solverOptions = osnspb.numericsSolverOptions()
                    fc_index=1
                    # Friction one-contact solver options
                    fc_internal_solver_options = solverOptions.internalSolvers[fc_index]
                    fc_internal_solver_options.iparam[0] = 100  # Local solver iterations
                    fc_internal_solver_options.solverId = local_solver
                else:
                    if 'NewtonImpactFrictionNSL' in set(nslaw_type_list):
                        if self._dimension ==3:
                            osnspb=FrictionContact(3, solver)
                            solverOptions = osnspb.numericsSolverOptions()
                            fc_index=0
                            # Friction one-contact solver options
                            fc_internal_solver_options = solverOptions.internalSolvers[fc_index]
                            fc_internal_solver_options.iparam[0] = 100  # Local solver iterations
                            fc_internal_solver_options.solverId = local_solver
                        elif self._dimension ==2:
                            osnspb=FrictionContact(2)
                        
                        osnspb.setMaxSize(osnspb_max_size)
                        osnspb.setMStorageType(1)
                    elif 'NewtonImpactRollingFrictionNSL' in set(nslaw_type_list):
                         osnspb=RollingFrictionContact(5)
                         solverOptions = osnspb.numericsSolverOptions()
                         solverOptions.iparam[0]=itermax
                         solverOptions.dparam[0] = tolerance
                    else:
                        raise RuntimeError("Unknown nslaw type"+ str(set(nslaw_type_list)))

        else:
            if (osi == Kernel.MoreauJeanGOSI):
                if (nb_of_nslaw_type >1) or 'NewtonImpactFrictionNSL' not in set(nslaw_type_list):
                    raise RuntimeError("MoreauJeanGOSI can deal only with NewtonImpactFrictionNSL nslaw")
                else:
                    from siconos.io.FrictionContactTrace import GlobalFrictionContactTrace
                    if (solver == Numerics.SICONOS_FRICTION_3D_NSGS or solver == Numerics.SICONOS_GLOBAL_FRICTION_3D_ADMM):
                        osnspb=GlobalFrictionContactTrace(3,
                                                          SICONOS_GLOBAL_FRICTION_3D_ADMM,
                                                          friction_contact_trace_params,nsds)
                        osnspb.setMStorageType(2)
                    else:
                        osnspb=GlobalFrictionContactTrace(3, solver,
                                                          friction_contact_trace_params,nsds)
                        solverOptions = osnspb.numericsSolverOptions()
                        # Friction one-contact solver options
                        fc_index=0
                        fc_internal_solver_options = solverOptions.internalSolvers[fc_index]
                        fc_internal_solver_options.iparam[0] = itermax
                        fc_internal_solver_options.dparam[0] = tolerance
                        osnspb.setMStorageType(2)
                    osnspb.setMaxSize(osnspb_max_size)
            else:
                from siconos.io.FrictionContactTrace import FrictionContactTrace
                osnspb=FrictionContactTrace(3, solver,friction_contact_trace_params,nsds)
                osnspb.setMaxSize(osnspb_max_size)
                osnspb.setMStorageType(1)


        # Numerics solver general  options
        solverOptions = osnspb.numericsSolverOptions()
        solverOptions.iparam[0] = itermax
        solverOptions.dparam[0] = tolerance

        # -- full error evaluation
        #solverOptions.iparam[1]=Numerics.SICONOS_FRICTION_3D_NSGS_ERROR_EVALUATION_FULL
        # --  Adaptive error evaluation
        #solverOptions.iparam[1]=Numerics.SICONOS_FRICTION_3D_NSGS_ERROR_EVALUATION_ADAPTIVE
        #solverOptions.iparam[8]=1
        # -- light error evaluation with full final
        if self._dimension ==3: # ugly
            solverOptions.iparam[1] = Numerics.SICONOS_FRICTION_3D_NSGS_ERROR_EVALUATION_LIGHT
            solverOptions.iparam[14] = Numerics.SICONOS_FRICTION_3D_NSGS_FILTER_LOCAL_SOLUTION_TRUE

        
        osnspb.setNumericsVerboseMode(numerics_verbose)
        if numerics_verbose:
            Numerics.numerics_set_verbose(numerics_verbose_level)

        # keep previous solution
        osnspb.setKeepLambdaAndYState(True)

        # (6) Simulation setup with (1) (2) (3) (4) (5)
        if time_stepping == Kernel.TimeSteppingDirectProjection:
            osnspb_pos=Kernel.MLCPProjectOnConstraints(Numerics.SICONOS_MLCP_ENUM, 1.0)
            so_pos = osnspb.numericsSolverOptions()
            so_pos.iparam[0]=itermax
            so_pos.dparam[0]=tolerance
            osnspb_pos.setMaxSize(osnspb_max_size)
            osnspb_pos.setMStorageType(0) # "not yet implemented for sparse storage"
            osnspb_pos.setNumericsVerboseMode(numerics_verbose)
            osnspb_pos.setKeepLambdaAndYState(True)
            simulation=time_stepping(nsds,timedisc, self._osi, osnspb, osnspb_pos)
            simulation.setProjectionMaxIteration(projection_itermax)
            simulation.setConstraintTolUnilateral(projection_tolerance_unilateral);
            simulation.setConstraintTol(projection_tolerance);
        else:
            simulation=time_stepping(nsds,timedisc)
            simulation.insertIntegrator(self._osi)

        simulation.insertNonSmoothProblem(osnspb)
        simulation.insertInteractionManager(self._interman)

        simulation.setNewtonOptions(Newton_options)
        simulation.setNewtonMaxIteration(Newton_max_iter)
        simulation.setNewtonTolerance(1e-10)
        if verbose:
            simulation.setDisplayNewtonConvergence(True)


        self._simulation = simulation
        # input()

        if len(self._plugins) > 0:
            self.print_verbose ('import plugins ...')
            self.import_plugins()



        if len(self._external_functions) > 0:
            self.print_verbose ('import external functions ...')
            self.import_external_functions()


        if controller is not None:
            controller.initialize(self)


        self.print_verbose ('first output static and dynamic objects ...')
        self.output_static_objects()
        self.output_dynamic_objects()

        if self._should_output_domains:
            self.log(self.output_domains, with_timer)()

        # nsds=model.nonSmoothDynamicalSystem()
        # nds= nsds.getNumberOfDS()
        # for i in range(nds):
        #     ds=nsds.dynamicalSystem(i)
        #     ds.display()
        # raw_input()
        self.print_verbose ('start simulation ...')
        self._initializing=False
        while simulation.hasNextEvent():

            if verbose_progress:
                self.print_verbose('step', k, 'of', k0 + int((T - t0) / h)-1)

            if start_run_iteration_hook is not None:
                self.log(start_run_iteration_hook, with_timer)(self)
                
            self.log(self.import_births, with_timer)(body_class,
                                                     shape_class,
                                                     face_class,
                                                     edge_class)

            self.log(self.execute_deaths, with_timer)()

            if controller is not None:
                controller.step()

            if (friction_contact_trace == True) :
                 osnspb._stepcounter = k


            if explode_Newton_solve:
                if(time_stepping == TimeStepping) :
                    self.log(self.explode_Newton_solve, with_timer, before=False)(with_timer)
                else:
                    print('simulation of type', type(time_stepping),' has no exploded version' )
                    self.log(simulation.computeOneStep, with_timer)()
            else:
                self.log(simulation.computeOneStep, with_timer)()


            if (self._output_frequency and (k % self._output_frequency == 0)) or (k == 1):
                if verbose:
                    self.print_verbose ('output in hdf5 file at step ', k)

                self.log(self.output_dynamic_objects, with_timer)()

                self.log(self.output_velocities, with_timer)()

                self.log(self.output_contact_forces, with_timer)()

                if self._should_output_domains:
                    self.log(self.output_domains, with_timer)()

                self.log(self.output_solver_infos, with_timer)()

                self.log(self._out.flush)()

            if self._output_backup:
                if (k % self._output_backup_frequency == 0) or (k == 1):
                    shutil.copyfile(self._io_filename, self._io_filename_backup)

            self.log(simulation.clearNSDSChangeLog, with_timer)()

            # Note these are not the same and neither is correct.
            # "_interman.statistics" gives the number of contacts
            # collected by the collision engine, but it's possible some
            # are not in indexset1.  Meanwhile checking the size of
            # the non-smooth problem is wrong when there are joints.
            if use_bullet:
                number_of_contacts = (
                    self._interman.statistics().new_interactions_created
                    + self._interman.statistics().existing_interactions_processed)
                if verbose and number_of_contacts > 0 :
                    self.print_verbose('number of contacts', number_of_contacts, '(detected)', osnspb.getSizeOutput()//3, '(active at velocity level. approx)')
                    self.print_solver_infos()

            else:
                number_of_contacts = osnspb.getSizeOutput()//3
                if verbose and number_of_contacts > 0 :
                    self.print_verbose('number of active contacts  at the velocity level (approx)', number_of_contacts)
                    self.print_solver_infos()
            

            if violation_verbose and number_of_contacts > 0 :
                if len(simulation.y(0,0)) >0 :
                    self.print_verbose('violation info')
                    y=simulation.y(0,0)
                    yplus=  np.zeros((2,len(y)))
                    yplus[0,:]=y
                    y=np.min(yplus,axis=1)
                    violation_max=np.max(-y)
                    self.print_verbose('  violation max :',violation_max)
                    if self._collision_margin is not None:
                        if  (violation_max >= self._collision_margin):
                            self.print_verbose('  violation max is larger than the collision_margin')
                    lam=simulation.lambda_(1,0)
                    self.print_verbose('  lambda max :',np.max(lam))
                    #print(' lambda : ',lam)
                    #raw_input()


                if len(simulation.y(1,0)) >0 :
                    v=simulation.y(1,0)
                    vplus=  np.zeros((2,len(v)))
                    vplus[0,:]=v
                    v=np.max(vplus,axis=1)
                    self.print_verbose('  velocity max :',np.max(v))
                    self.print_verbose('  velocity min :',np.min(v))
                #     #print(simulation.output(1,0))

            precision = solverOptions.dparam[Numerics.SICONOS_DPARAM_RESIDU]
            if (exit_tolerance is not None):
                if (precision > exit_tolerance):
                    print('precision is larger exit_tolerance')
                    return False
            self.log(simulation.nextStep, with_timer)()

            if end_run_iteration_hook is not None:
                self.log(end_run_iteration_hook, with_timer) (self)

            self.print_verbose ('')
            k += 1
        return True
