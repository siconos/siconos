# Mechanics IO

from __future__ import print_function

import os
import sys

from math import cos, sin

import shlex
import numpy as np
import h5py
import bisect

import tempfile
from contextlib import contextmanager
from operator import itemgetter

try:
    import vtk
except:
    pass

from siconos.mechanics.contact_detection.tools import Contactor

from siconos.mechanics import joints

try:
    from siconos.mechanics.contact_detection.bullet import \
        BulletDS, BulletWeightedShape, \
        btCollisionObject, btQuaternion, btTransform, btVector3, quatRotate

    from siconos.mechanics.contact_detection.bullet import \
        cast_BulletR

    from siconos.mechanics.contact_detection.bullet import \
        __mul__ as mul

    from siconos.mechanics.contact_detection.bullet import btVector3, \
        btConvexHullShape, btCylinderShape, btBoxShape, btSphereShape, \
        btConeShape, btCapsuleShape, btCompoundShape, btTriangleIndexVertexArray, \
        btGImpactMeshShape

    bullet_is_here = True

except:

    bullet_is_here = False

    pass


try:
    from siconos.mechanics.occ import \
        OccContactShape, OccBody, OccContactFace, OccContactEdge, \
        OccTimeStepping, OccSpaceFilter, OccWrap
except:
    pass

from siconos.kernel import \
    cast_NewtonImpactFrictionNSL, EqualityConditionNSL, Interaction

import siconos.kernel as Kernel
from siconos.io.io_base import MechanicsIO

import siconos.numerics as Numerics

from scipy import constants

import time


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


def group(h, name):
    try:
        return h[name]
    except KeyError:
        return h.create_group(name)


def data(h, name, nbcolumns):
    try:
        return h[name]
    except KeyError:
        return h.create_dataset(name, (0, nbcolumns),
                                maxshape=(None, nbcolumns))


def add_line(dataset, line):
    dataset.resize(dataset.shape[0] + 1, 0)
    dataset[dataset.shape[0] - 1, :] = line


def str_of_file(filename):
    with open(filename, 'r') as f:
        return str(f.read())


#
# load .vtp file
#
def loadMesh(shape_filename, collision_margin):
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
        shape.setMargin(collision_margin)
        for p in coors:
                shape.addPoint(btVector3(*p))

    return keep, shape


class ShapeCollection():

    """
    Instantiation of added contact shapes
    """

    def __init__(self, io, collision_margin=0.04):
        self._io = io
        self._shapes = dict()
        self._tri = dict()
        self._collision_margin = collision_margin
        # print('self._collision_margin',self._collision_margin)
        if bullet_is_here:

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
            edge_class=None):

        if not shape_name in self._shapes:

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
                            (self._tri[shape_name],
                             self._shapes[shape_name]) = loadMesh(
                                 tmpf[1], self._collision_margin)
                    else:
                        assert False
                elif self.attributes(shape_name)['type'] in['step']:
                    from OCC.STEPControl import STEPControl_Reader
                    from OCC.BRep import BRep_Builder
                    from OCC.TopoDS import TopoDS_Compound
                    from OCC.IFSelect import IFSelect_RetDone,\
                        IFSelect_ItemsByEntity

                    builder = BRep_Builder()
                    comp = TopoDS_Compound()
                    builder.MakeCompound(comp)

                    step_name = self.attributes(shape_name)['step']

                    with tmpfile(contents=self.shape(step_name)[:][0]) as tmpf:
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

                            l = []

                            for i in range(1, nbs + 1):
                                shape = step_reader.Shape(i)
                                builder.Add(comp, shape)

                            self._shapes[shape_name] = comp
                            self._io._keep.append(self._shapes[shape_name])

                elif self.attributes(shape_name)['type'] in['brep']:
                    if not 'contact' in self.attributes(shape_name):

                        # the reference brep
                        if shape_class is None:
                            brep_class = OccContactShape
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
                        assert 'index' in self.attributes(shape_name)
                        assert 'brep' in self.attributes(shape_name)
                        contact_index = self.attributes(shape_name)['index']

                        ref_brep = self.get(
                            self.attributes(shape_name)['brep'], shape_class)

                        if self.attributes(shape_name)['contact'] == 'Face':
                            if face_class is None:
                                face_maker = OccContactFace
                            else:
                                face_maker = face_class

                            self._shapes[shape_name] = \
                                face_maker(ref_brep,
                                           contact_index)

                        elif self.attributes(shape_name)['contact'] == 'Edge':
                            if edge_class is None:
                                edge_maker = OccContactEdge
                            else:
                                edge_maker = edge_class
                            self._shapes[shape_name] = \
                                edge_maker(ref_brep,
                                           contact_index)

                        self._io._keep.append(self._shapes[shape_name])
                else:
                    # a convex point set
                    convex = btConvexHullShape()
                    convex.setMargin(self._collision_margin)
                    for points in self.shape(shape_name):
                        convex.addPoint(btVector3(float(points[0]),
                                                  float(points[1]),
                                                  float(points[2])))
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
                    self._shapes[shape_name] = primitive(
                        btVector3(attrs[0] / 2,
                                  attrs[
                                  1] / 2,
                                  attrs[2] / 2))
                elif name in ['Cylinder']:
                    self._shapes[shape_name] = primitive(btVector3(attrs[0],
                                                                   attrs[
                                                                       1] / 2,
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
                    self._shapes[shape_name] = primitive(*attrs)

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
                 broadphase=None, osi=None, shape_filename=None,
                 set_external_forces=None, gravity_scale=None, collision_margin=None):

        if io_filename is None:
            self._io_filename = '{0}.hdf5'.format(
                os.path.splitext(os.path.basename(sys.argv[0]))[0])
        else:
            self._io_filename = io_filename
        self._mode = mode
        self._broadphase = broadphase
        self._osi = osi
        self._static_origins = []
        self._static_orientations = []
        self._static_transforms = []
        self._static_cobjs = []
        self._shape = None
        self._shapeid = dict()
        self._static_data = None
        self._velocities_data = None
        self._dynamic_data = None
        self._cf_data = None
        self._solv_data = None
        self._input = None
        self._nslaws = None
        self._out = None
        self._data = None
        self._joints = None
        self._io = MechanicsIO()
        self._set_external_forces = set_external_forces
        self._shape_filename = shape_filename
        self._number_of_shapes = 0
        self._number_of_dynamic_objects = 0
        self._number_of_static_objects = 0
        self._gravity_scale = gravity_scale
        self._collision_margin = collision_margin
        self._output_frequency = 1
        self._keep = []
        self._scheduled_births = []
        self._births = dict()

    def __enter__(self):
        if self._set_external_forces is None:
            self._set_external_forces = self.apply_gravity

        if self._gravity_scale is None:
            self._gravity_scale = 1  # 1 => m, 1/100. => cm

        self._out = h5py.File(self._io_filename, self._mode)
        self._data = group(self._out, 'data')
        self._ref = group(self._data, 'ref')
        self._joints = group(self._data, 'joints')
        self._static_data = data(self._data, 'static', 9)
        self._velocities_data = data(self._data, 'velocities', 8)
        self._dynamic_data = data(self._data, 'dynamic', 9)
        self._cf_data = data(self._data, 'cf', 15)
        self._solv_data = data(self._data, 'solv', 4)
        self._input = group(self._data, 'input')
        self._nslaws = group(self._data, 'nslaws')

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
        return self._nslaws

    def joints(self):
        """
        Joints between dynamic objects or between an object and the scenery.
        """
        return self._joints

    def importNonSmoothLaw(self, name):
        if self._broadphase is not None:
            nslawClass = getattr(Kernel, self._nslaws[name].attrs['type'])
            # only this one at the moment
            assert(nslawClass == Kernel.NewtonImpactFrictionNSL)
            nslaw = nslawClass(float(self._nslaws[name].attrs['e']), 0.,
                               float(self._nslaws[name].attrs['mu']), 3)
            self._broadphase.insert(nslaw,
                                    int(self._nslaws[name].attrs['gid1']),
                                    int(self._nslaws[name].attrs['gid2']))

    def importBRepObject(self, name, translation, orientation,
                         velocity, contactors, mass, given_inertia, body_class,
                         shape_class, face_class, edge_class):

        if body_class is None:
            body_class = OccBody

        if given_inertia is not None:
            inertia = given_inertia
        else:
            inertia = np.eye(3)

        if mass == 0.:
            # a static object
            pass

        else:
            body = body_class(
                translation + orientation, velocity, mass, inertia)
            for contactor in contactors:

                # /!\ shared pointer <-> python ref ...
                shape_instantiated = self._shape.get(contactor.name,
                                                     shape_class, face_class, edge_class)
                self._keep.append(shape_instantiated)

                body.addContactShape(shape_instantiated,
                                     contactor.translation,
                                     contactor.orientation,
                                     contactor.group)
            self._set_external_forces(body)

            # add the dynamical system to the non smooth
            # dynamical system
            nsds = self._broadphase.model().nonSmoothDynamicalSystem()
            nsds.topology().setOSI(body, self._osi)
            nsds.setName(body, str(name))

    def importBulletObject(self, name, translation, orientation,
                           velocity, contactors, mass, inertia,
                           body_class, shape_class, birth=False):

        if body_class is None:
            body_class = BulletDS

        if self._broadphase is not None and 'input' in self._data:
            if mass == 0.:
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

                    self._static_origins.append(
                        static_cobj.getWorldTransform().getOrigin())

                    self._static_orientations.append(
                        static_cobj.getWorldTransform().getRotation())

                    self._static_transforms.append(
                        static_cobj.getWorldTransform())

                    static_cobj.setCollisionShape(
                        self._shape.get(c.name))
                    self._static_cobjs.append(static_cobj)

                    self._broadphase.addStaticObject(static_cobj,
                                                     int(c.group))

            else:
                # a moving object
                bws = BulletWeightedShape(
                    self._shape.get(contactors[0].name), mass)

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

                for contactor in contactors[1:]:
                    shape_id = self._shapeid[contactor.name]

                    body.addCollisionShape(self._shape.get(contactor.name),
                                           contactor.translation,
                                           contactor.orientation,
                                           contactor.group)

                # set external forces
                self._set_external_forces(body)

                # add the dynamical system to the non smooth
                # dynamical system
                if birth:
                    self._broadphase.addDynamicObject(
                        body,
                        self._broadphase.model().simulation(),
                        self._osi)
                    nsds = self._broadphase.model().nonSmoothDynamicalSystem()
                    nsds.setName(body, str(name))
                else:
                    nsds = self._broadphase.model().nonSmoothDynamicalSystem()
                    nsds.insertDynamicalSystem(body)
                    nsds.topology().setOSI(body, self._osi)
                    nsds.setName(body, str(name))

    def importJoint(self, name):
        if self._broadphase is not None:
            topo = self._broadphase.model().nonSmoothDynamicalSystem().\
                topology()

            joint_class = getattr(joints,
                                  self.joints()[name].attrs['type'])

            joint_nslaw = EqualityConditionNSL(5)

            ds1_name = self.joints()[name].attrs['object1']
            ds1 = topo.getDynamicalSystem(ds1_name)

            if 'object2' in self.joints()[name].attrs:
                ds2_name = self.joints()[name].attrs['object2']
                ds2 = topo.getDynamicalSystem(ds2_name)
                joint = joint_class(ds1,
                                    ds2,
                                    self.joints()[name].attrs['pivot_point'],
                                    self.joints()[name].attrs['axis'])
                joint_inter = Interaction(5, joint_nslaw, joint)
                self._broadphase.model().nonSmoothDynamicalSystem().\
                    link(joint_inter, ds1, ds2)

            else:
                joint = joint_class(ds1,
                                    self.joints()[name].attrs['pivot_point'],
                                    self.joints()[name].attrs['axis'])

                joint_inter = Interaction(5, joint_nslaw, joint)
                self._broadphase.model().nonSmoothDynamicalSystem().\
                    link(joint_inter, ds1)

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

        for shape_name in self._ref:
            self._shapeid[shape_name] = self._ref[shape_name].attrs['id']
            self._number_of_shapes += 1

        # import dynamical systems
        if self._broadphase is not None and 'input' in self._data:
            counter =0
            for (name, obj) in sorted(self._input.items(),
                                      key=lambda x: x[0]):
                counter +=1
                print ('import  object ', counter)
                input_ctrs = [ctr for _n_, ctr in obj.items()]
                mass = obj.attrs['mass']
                time_of_birth = obj.attrs['time_of_birth']

                if time_of_birth > time:
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
                    if mass > 0 and self.dynamic_data() is not None and \
                            len(self.dynamic_data()) > 0:
                        dpos_data = self.dynamic_data()
                        max_time = max(dpos_data[:, 0])
                        id_last = np.where(
                            abs(dpos_data[:, 0] - max_time) < 1e-9)[0]
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
                        velocity = (xvel[2], xvel[3], xvel[4])

                    # start from initial conditions
                    else:
                        translation = obj.attrs['translation']
                        orientation = obj.attrs['orientation']
                        velocity = obj.attrs['velocity']

                    contactors = [Contactor(ctr.attrs['name'],
                                            int(ctr.attrs['group']),
                                            floatv(ctr.attrs['translation']),
                                            floatv(ctr.attrs['orientation']))
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
                        self.importBRepObject(
                            name, floatv(translation), floatv(orientation),
                            floatv(
                                velocity), contactors, float(mass),
                            inertia, body_class, shape_class, face_class,
                            edge_class)
                    else:
                        # Bullet object
                        self.importBulletObject(
                            name, floatv(translation), floatv(orientation),
                            floatv(velocity), contactors, float(mass),
                            inertia, body_class, shape_class)

            # import nslaws
            # note: no time of birth for nslaws and joints
            for name in self._nslaws:
                self.importNonSmoothLaw(name)

            for name in self.joints():
                self.importJoint(name)

    def importBirths(self, body_class=None, shape_class=None,
                     face_class=None, edge_class=None,):
        """
        Import new objects in the broadphase.
        """
        time = self._broadphase.model().simulation().nextTime()

        ind_time = bisect.bisect_left(self._scheduled_births, time)

        current_times_of_births = self._scheduled_births[:ind_time]
        self._scheduled_births = self._scheduled_births[ind_time:]

        #print (time, current_times_of_births)
        for time_of_birth in current_times_of_births:

            for (name, obj) in self._births[time_of_birth]:
                translation = obj.attrs['translation']
                orientation = obj.attrs['orientation']
                velocity = obj.attrs['velocity']

                input_ctrs = [ctr for _n_, ctr in obj.items()]
                mass = obj.attrs['mass']

                contactors = [Contactor(ctr.attrs['name'],
                                        int(ctr.attrs['group']),
                                        floatv(ctr.attrs['translation']),
                                        floatv(ctr.attrs['orientation']))
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
                    self.importBRepObject(
                        name, floatv(translation), floatv(orientation),
                        floatv(
                            velocity), contactors, float(mass),
                        inertia, body_class, shape_class, face_class,
                        edge_class)
                else:
                    # Bullet object
                    self.importBulletObject(
                        name, floatv(translation), floatv(orientation),
                        floatv(velocity), contactors, float(mass),
                        inertia, body_class, shape_class, birth=True)

    def outputStaticObjects(self):
        """
        Outputs translations and orientations of static objects
        """
        time = self._broadphase.model().simulation().nextTime()
        idd = -1
        p = 0
        self._static_data.resize(len(self._static_transforms), 0)

        for transform in self._static_transforms:
            translation = transform.getOrigin()
            rotation = transform.getRotation()
            self._static_data[p, :] = \
                [time,
                 idd,
                 translation.x(),
                 translation.y(),
                 translation.z(),
                 rotation.w(),
                 rotation.x(),
                 rotation.y(),
                 rotation.z()]
            idd -= 1
            p += 1

    def outputDynamicObjects(self):
        """
        Outputs translations and orientations of dynamic objects.
        """

        current_line = self._dynamic_data.shape[0]

        time = self._broadphase.model().simulation().nextTime()

        positions = self._io.positions(self._broadphase.model())

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

        time = self._broadphase.model().simulation().nextTime()

        velocities = self._io.velocities(self._broadphase.model())

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
        """
        if self._broadphase.model().nonSmoothDynamicalSystem().\
                topology().indexSetsSize() > 1:
            time = self._broadphase.model().simulation().nextTime()
            contact_points = self._io.contactPoints(self._broadphase.model())

            if contact_points is not None:

                current_line = self._cf_data.shape[0]
                self._cf_data.resize(current_line + contact_points.shape[0], 0)
                times = np.empty((contact_points.shape[0], 1))
                times.fill(time)

                self._cf_data[current_line:, :] = \
                    np.concatenate((times,
                                    contact_points),
                                   axis=1)

    def outputSolverInfos(self):
        """
        Outputs solver #iterations & precision reached
        """

        time = self._broadphase.model().simulation().nextTime()
        so = self._broadphase.model().simulation().oneStepNSProblem(0).\
            numericsSolverOptions()

        current_line = self._solv_data.shape[0]
        self._solv_data.resize(current_line + 1, 0)
        if so.solverId == Numerics.SICONOS_GENERIC_MECHANICAL_NSGS:
            iterations = so.iparam[3]
            precision = so.dparam[2]
            local_precision = so.dparam[3]
        elif so.solverId == Numerics.SICONOS_FRICTION_3D_NSGS:
            iterations = so.iparam[7]
            precision = so.dparam[1]
            local_precision = 0.
        # maybe wrong for others
        else:
            iterations = so.iparam[1]
            precision = so.dparam[1]
            local_precision = so.dparam[2]

        self._solv_data[current_line, :] = [time, iterations, precision,
                                            local_precision]

    def printSolverInfos(self):
        """
        Outputs solver #iterations & precision reached
        """
        time = self._broadphase.model().simulation().nextTime()
        so = self._broadphase.model().simulation().oneStepNSProblem(0).\
            numericsSolverOptions()
        if so.solverId == Numerics.SICONOS_GENERIC_MECHANICAL_NSGS:
            iterations = so.iparam[3]
            precision = so.dparam[2]
            local_precision = so.dparam[3]
        elif so.solverId == Numerics.SICONOS_FRICTION_3D_NSGS:
            iterations = so.iparam[7]
            precision = so.dparam[1]
            local_precision = 0.
        # maybe wrong for others
        else:
            iterations = so.iparam[1]
            precision = so.dparam[1]
            local_precision = so.dparam[2]

        
        print('SolverInfos at time :', time,
              'iterations= ', iterations,
              'precision=', precision,
              'local_precision=', )

    def addMeshFromString(self, name, shape_data):
        """
        Add a mesh shape from a string.
        Accepted format : mesh encoded in VTK .vtp format
        """
        if name not in self._ref:

            shape = self._ref.create_dataset(name, (1,),
                                             dtype=h5py.new_vlen(str))
            shape[:] = shape_data
            shape.attrs['id'] = self._number_of_shapes
            shape.attrs['type'] = 'vtp'
            self._shapeid[name] = shape.attrs['id']
            self._number_of_shapes += 1

    def addMeshFromFile(self, name, filename):
        """
        Add a mesh shape from a file.
        Accepted format : .stl or mesh encoded in VTK .vtp format
        """
        if name not in self._ref:

            if os.path.splitext(filename)[-1][1:] == 'stl':
                reader = vtk.vtkSTLReader()
                reader.SetFileName(filename)
                reader.Update()

                with tmpfile() as tmpf:
                    writer = vtk.vtkXMLPolyDataWriter()
                    writer.SetInputData(reader.GetOutput())
                    writer.SetFileName(tmpf[1])
                    writer.Write()

                    shape_data = str_of_file(tmpf[1])

            else:
                assert os.path.splitext(filename)[-1][1:] == 'vtp'
                shape_data = str_of_file(filename)

            self.addMeshFromString(name, shape_data)

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
        Add an OpenCascade TopoDS_Shape
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

    def addContactFromBRep(self, name, brepname, contact_type,
                           index, collision_group=0, associated_shape=None):
        """
        Add contact reference from a previously added brep.
        """
        if name not in self._ref:
            shape = self._ref.create_dataset(name, (1,),
                                             dtype=h5py.new_vlen(str))
            shape.attrs['id'] = self._number_of_shapes
            shape.attrs['type'] = 'brep'
            shape.attrs['contact'] = contact_type
            shape.attrs['brep'] = brepname
            shape.attrs['index'] = index
            if associated_shape is not None:
                shape.attrs['associated_shape'] = associated_shape
            self._shapeid[name] = shape.attrs['id']
            self._number_of_shapes += 1

    def addContactFromOccShape(self, name, occ_shape_name, contact_type,
                               index, collision_group=0, associated_shape=None):
        """
        Add contact reference from a previously added brep.
        """
        if name not in self._ref:
            shape = self._ref.create_dataset(name, (1,),
                                             dtype=h5py.new_vlen(str))
            shape.attrs['id'] = self._number_of_shapes
            shape.attrs['type'] = 'step'
            shape.attrs['contact'] = contact_type
            shape.attrs['step'] = occ_shape_name
            shape.attrs['index'] = index
            if associated_shape is not None:
                shape.attrs['associated_shape'] = associated_shape
            self._shapeid[name] = shape.attrs['id']
            self._number_of_shapes += 1

    def addConvexShape(self, name, points):
        """
        Add a convex shape defined by a list of points.
        """
        if name not in self._ref:
            apoints = np.array(points)
            shape = self._ref.create_dataset(name,
                                             (apoints.shape[0],
                                              apoints.shape[1]))
            shape[:] = points[:]
            shape.attrs['type'] = 'convex'
            shape.attrs['id'] = self._number_of_shapes
            self._shapeid[name] = shape.attrs['id']
            self._number_of_shapes += 1

    def addPrimitiveShape(self, name, primitive, params):
        """
        Add a primitive shape.
        """
        if name not in self._ref:
            shape = self._ref.create_dataset(name, (1, len(params)))
            shape.attrs['id'] = self._number_of_shapes
            shape.attrs['type'] = 'primitive'
            shape.attrs['primitive'] = primitive
            shape[:] = params
            self._shapeid[name] = shape.attrs['id']
            self._number_of_shapes += 1

    def addObject(self, name, shapes,
                  translation,
                  orientation=[1, 0, 0, 0],
                  velocity=[0, 0, 0, 0, 0, 0],
                  mass=0, inertia=None, time_of_birth=0.):
        """
        Add an object with associated contact shapes that are called
        contactors.
        Contact detection is defined by the list of contactors.
        The initial translation is mandatory : [x, y z].
        If the mass is zero this is a static object.
        A default inertia is associated if it is not given.
        A time of birth may be specified, so the object will not appear
        in the scene until this time.
        """
        # print(arguments())
        if len(orientation) == 2:
            # axis + angle
            axis = orientation[0]
            assert len(axis) == 3
            angle = orientation[1]
            assert type(angle) is float
            n = sin(angle / 2.) / np.linalg.norm(axis)

            ori = [cos(angle / 2.), axis[0] * n, axis[1] * n, axis[2] * n]
        else:
            # a given quaternion
            ori = orientation

        if name not in self._input:

            obj = group(self._input, name)

            obj.attrs['time_of_birth'] = time_of_birth

            obj.attrs['mass'] = mass
            obj.attrs['translation'] = translation
            obj.attrs['orientation'] = ori
            obj.attrs['velocity'] = velocity

            if inertia is not None:
                obj.attrs['inertia'] = inertia

            for num, shape in enumerate(shapes):
                dat = data(obj, '{0}-{1}'.format(shape.name, num), 0)
                dat.attrs['name'] = shape.name
                if hasattr(shape, 'group'):
                    dat.attrs['group'] = shape.group
                if hasattr(shape, 'parameters') and \
                        shape.parameters is not None:
                    dat.attrs['parameters'] = shape.parameters
                dat.attrs['translation'] = shape.translation
                dat.attrs['orientation'] = shape.orientation

            if mass == 0:
                obj.attrs['id'] = - (self._number_of_static_objects + 1)
                self._number_of_static_objects += 1

            else:
                obj.attrs['id'] = (self._number_of_dynamic_objects + 1)
                self._number_of_dynamic_objects += 1

    def addNewtonImpactFrictionNSL(self, name, mu, e=0, collision_group1=0,
                                   collision_group2=0):
        """
        Add a nonsmooth law for contact between 2 groups.
        Only NewtonImpactFrictionNSL are supported.
        name is a user identifiant and must be unique,
        mu is the coefficient of friction,
        e is the coefficient of restitution on the contact normal,
        gid1 and gid2 define the group identifiants.

        """
        if name not in self._nslaws:
            nslaw = self._nslaws.create_dataset(name, (0,))
            nslaw.attrs['type'] = 'NewtonImpactFrictionNSL'
            nslaw.attrs['mu'] = mu
            nslaw.attrs['e'] = e
            nslaw.attrs['gid1'] = collision_group1
            nslaw.attrs['gid2'] = collision_group2

    def addJoint(self, name, object1, object2=None, pivot_point=[0, 0, 0],
                 axis=[0, 1, 0],
                 joint_class='PivotJointR'):
        """
        add a pivot joint between two objects
        """
        if name not in self.joints():
            joint = self.joints().create_dataset(name, (0,))
            joint.attrs['object1'] = object1
            if object2 is not None:
                joint.attrs['object2'] = object2
            joint.attrs['type'] = joint_class
            joint.attrs['pivot_point'] = pivot_point
            joint.attrs['axis'] = axis

    def run(self,
            with_timer=False,
            time_stepping=None,
            space_filter=None,
            body_class=None,
            shape_class=None,
            face_class=None,
            edge_class=None,
            gravity_scale=1.0,
            t0=0,
            T=10,
            h=0.0005,
            multipoints_iterations=True,
            theta=0.50001,
            Newton_max_iter=20,
            set_external_forces=None,
            solver=Numerics.SICONOS_FRICTION_3D_NSGS,
            itermax=100000,
            tolerance=1e-8,
            numerics_verbose=False,
            violation_verbose=False,
            output_frequency=None):
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

        """
        print ('load siconos module ...')
        from siconos.kernel import \
            Model, NonSmoothDynamicalSystem, OneStepNSProblem, MoreauJeanOSI,\
            TimeDiscretisation, GenericMechanical, FrictionContact,\
            NewtonImpactFrictionNSL

        from siconos.numerics import SICONOS_FRICTION_3D_ONECONTACT_NSN_AC

        from siconos.mechanics.contact_detection.bullet import \
            btConvexHullShape, btCollisionObject, \
            btBoxShape, btQuaternion, btTransform, btConeShape, \
            BulletSpaceFilter, cast_BulletR, \
            BulletWeightedShape, BulletDS, BulletTimeStepping
        
        print ('setup model simulation ...')
        if set_external_forces is not None:
            self._set_external_forces = set_external_forces

        if time_stepping is None:
            time_stepping = BulletTimeStepping

        if space_filter is None:
            space_filter = BulletSpaceFilter

        if output_frequency is not None:
            self._output_frequency = output_frequency

        if gravity_scale is not None:
            self._gravity_scale = gravity_scale

        # cold restart
        times = set()
        if self.dynamic_data() is not None and len(self.dynamic_data()) > 0:
            dpos_data = self.dynamic_data()
            times = set(dpos_data[:, 0])
            t0 = float(max(times))
            T = float(t0 + T)
            print ('restart from previous simulation at t0={0}'.format(t0))
            print ('run until T={0}'.format(T))

        # Model
        #
        model = Model(t0, T)

        # (1) OneStepIntegrators
        joints = list(self.joints())

        self._osi = MoreauJeanOSI(theta)

        # (2) Time discretisation --
        timedisc = TimeDiscretisation(t0, h)

        if len(joints) > 0:
            osnspb = GenericMechanical(SICONOS_FRICTION_3D_ONECONTACT_NSN_AC)
        else:
            osnspb = FrictionContact(3, solver)

        osnspb.numericsSolverOptions().iparam[0] = itermax
        osnspb.numericsSolverOptions().internalSolvers.iparam[0] = 100
        osnspb.numericsSolverOptions().dparam[0] = tolerance
        osnspb.setMaxSize(30000)
        osnspb.setMStorageType(1)
        osnspb.setNumericsVerboseMode(numerics_verbose)

        # keep previous solution
        osnspb.setKeepLambdaAndYState(True)

        # (5) broadphase contact detection
        self._broadphase = space_filter(model)
        if not multipoints_iterations:
            print("""
            ConvexConvexMultipointIterations and PlaneConvexMultipointIterations are unset
            """)
        else:
            if hasattr(self._broadphase, 'collisionConfiguration'):
                self._broadphase.collisionConfiguration().\
                    setConvexConvexMultipointIterations()
                self._broadphase.collisionConfiguration().\
                    setPlaneConvexMultipointIterations()

        # (6) Simulation setup with (1) (2) (3) (4) (5)
        simulation = time_stepping(timedisc)
        simulation.insertIntegrator(self._osi)
        simulation.insertNonSmoothProblem(osnspb)
        simulation.setNewtonMaxIteration(Newton_max_iter)
        simulation.setNewtonTolerance(1e-10)

        k0 = 1 + len(times)
        k = k0
        print ('import scene ...')
        self.importScene(t0, body_class, shape_class, face_class, edge_class)

        model.setSimulation(simulation)
        model.initialize()
        print ('first output static and dynamic objects ...')
        self.outputStaticObjects()
        self.outputDynamicObjects()

        # nsds = model.nonSmoothDynamicalSystem()
        # nds= nsds.getNumberOfDS()
        # for i in range(nds):
        #     ds = nsds.dynamicalSystem(i)
        #     ds.display()
        # raw_input()
        print ('start simulation ...')
        while simulation.hasNextEvent():

            print ('step', k, '<', k0 - 1 + int((T - t0) / h))

            log(self.importBirths(body_class=body_class,
                                  shape_class=shape_class,
                                  face_class=face_class,
                                  edge_class=edge_class))

            log(self._broadphase.buildInteractions, with_timer)\
                (model.currentTime())

            log(simulation.computeOneStep, with_timer)()

            if (k % self._output_frequency == 0) or (k == 1):
                print ('output in hdf5 file at step ', k)

                log(self.outputDynamicObjects, with_timer)()

                log(self.outputVelocities, with_timer)()

                log(self.outputContactForces, with_timer)()

                log(self.outputSolverInfos, with_timer)()

                log(self._out.flush)()


            numberOfContact=self._broadphase.model().simulation().oneStepNSProblem(0).getSizeOutput()/3
            if numberOfContact > 0 :
                print('number of contact',self._broadphase.model().simulation().oneStepNSProblem(0).getSizeOutput()/3)
                self.printSolverInfos()
                
            if violation_verbose and numberOfContact > 0 :
                if len(simulation.y(0,0)) >0 :
                    print('violation info')
                    y = simulation.y(0,0)
                    yplus=  np.zeros((2,len(y)))
                    yplus[0,:] = y
                    y=np.min(yplus,axis=1)
                    violation_max=np.max(-y)
                    print('  violation max :',violation_max)
                    if  (violation_max >= self._collision_margin):
                        print('  violation max is larger than the collision_margin')
                    lam = simulation.lambda_(1,0)
                    print('  lambda max :',np.max(lam))
                    #print(' lambda : ',lam)
                    #raw_input()


                if len(simulation.y(1,0)) >0 :
                    v = simulation.y(1,0)
                    vplus=  np.zeros((2,len(v)))
                    vplus[0,:] = v
                    v=np.max(vplus,axis=1)
                    print('  velocity max :',np.max(v))
                    print('  velocity min :',np.min(v))
                #     #print(simulation.output(1,0))


            log(simulation.nextStep, with_timer)()

            print ('')
            k += 1
