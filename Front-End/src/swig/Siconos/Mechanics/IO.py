# Mechanics IO

from __future__ import print_function

import os
import sys

from math import cos, sin

import shlex
import numpy as np
import h5py

import tempfile
from contextlib import contextmanager
from operator import itemgetter

try:
    import vtk
except:
    pass

from Siconos.Mechanics.ContactDetection import Contactor

from Siconos.Mechanics import Joints

try:
    from Siconos.Mechanics.ContactDetection.Bullet import \
        BulletDS, BulletWeightedShape, \
        btCollisionObject, btQuaternion, btTransform, btVector3, quatRotate

    from Siconos.Mechanics.ContactDetection.Bullet import \
        cast_BulletR

    from Siconos.Mechanics.ContactDetection.Bullet.BulletWrap import \
        __mul__ as mul

    from ContactDetection.Bullet import btVector3, \
        btConvexHullShape, btCylinderShape, btBoxShape, btSphereShape, \
        btConeShape, btCapsuleShape, btCompoundShape, btTriangleIndexVertexArray, \
        btGImpactMeshShape

except:
    pass


try:
    from Siconos.Mechanics.Occ import \
        OccContactShape, OccBody, OccContactFace, OccContactEdge, \
        OccTimeStepping, OccSpaceFilter
except:
    pass

from Siconos.Kernel import \
    cast_NewtonImpactFrictionNSL, EqualityConditionNSL, Interaction

import Siconos.Kernel as Kernel

from Siconos.IO import MechanicsIO

import Siconos.Numerics as Numerics

from scipy import constants

import time


@contextmanager
def tmpfile(suffix='', prefix='siconos_io'):
    """
    A context manager for a named temporary file.
    """
    (_, tfilename) = tempfile.mkstemp(suffix=suffix, prefix=prefix)
    fid = open(tfilename, 'w')
    yield (fid, tfilename)
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
    weight = [0, 0, - body.massValue() * g]
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
    dataset[dataset.shape[0]-1, :] = line


def str_of_file(filename):
    with open(filename, 'r') as f:
        return str(f.read())


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


class ShapeCollection():

    """
    Instantiation of added contact shapes
    """

    def __init__(self, io):
        self._io = io
        self._shapes = dict()

        self._primitive = {'Cylinder': btCylinderShape,
                           'Sphere': btSphereShape,
                           'Box': btBoxShape,
                           'Cone': btConeShape,
                           'Compound': btCompoundShape,
                           'Capsule': btCapsuleShape}

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


    def get(self, shape_name, shape_class=None, face_class=None, edge_class=None):

        if not shape_name in self._shapes:

            # load shape if it is an existing file
            if not isinstance(self.url(shape_name), str) and \
               not 'primitive' in self.attributes(shape_name):
                # assume a vtp file (xml) stored in a string buffer
                if self.attributes(shape_name)['type'] == 'vtk':
                    if self.shape(shape_name).dtype == h5py.new_vlen(str):
                        with tmpfile() as tmpf:
                            data = self.shape(shape_name)[:][0]
                            tmpf[0].write(data)
                            tmpf[0].flush()
                            self._tri[index], self._shape[index] = loadMesh(
                                tmpf[1])
                    else:
                        assert False
                elif self.attributes(shape_name)['type'] == 'brep':
                    if not 'contact' in self.attributes(shape_name):
                        # the reference brep
                        if shape_class is None:
                            brep = OccContactShape()
                        else:
                            brep = shape_class()

                        brep.importBRepFromString(self.shape(shape_name)[:][0])
                        self._shapes[shape_name] = brep
                        self._io._keep.append(self._shapes[shape_name])

                    else:
                        # a contact on a brep
                        assert 'contact' in self.attributes(shape_name)
                        assert 'index' in self.attributes(shape_name)
                        assert 'brep' in self.attributes(shape_name)
                        contact_index = self.attributes(shape_name)['index']
                        print (self.attributes(shape_name)['brep'])
                        ref_brep = self.get(self.attributes(shape_name)['brep'], shape_class)
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
                    for points in self.shape(shape_name):
                        convex.addPoint(btVector3(float(points[0]),
                                                  float(points[1]),
                                                  float(points[2])))
                    self._shapes[shape_name] = convex

            elif isinstance(self.url(shape_name), str) and \
                os.path.exists(self.url(shape_name)):
                self._tri[shape_name], self._shapes[shape_name] = loadMesh(
                    self.url(shape_name))
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
                    self._shapes[shape_name] = primitive(btVector3(attrs[0] / 2,
                                                                   attrs[1] / 2,
                                                                   attrs[2] / 2))
                elif name in ['Cylinder']:
                    self._shapes[shape_name] = primitive(btVector3(attrs[0],
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
                    self._shapes[shape_name] = primitive(*attrs)

        return self._shapes[shape_name]


class Hdf5():
    """a Hdf5 context manager reads at instantiation the positions and
       orientations of collision objects from hdf5 file

       It provides functions to output positions and orientations in
       the same file during simulation (output is done by default in
       pos.dat)

       with:
         time : float
         object_id : the object id (int)
         px, py, pz : components of the position (float)
         ow, ox, oy oz : components of an unit quaternion (float)

    """

    def __init__(self, io_filename = None, mode = 'w',
                 broadphase=None, osi=None, shape_filename=None,
                 set_external_forces=None, length_scale=None):

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
        self._length_scale = length_scale
        self._keep = []

    def __enter__(self):
        if self._set_external_forces is None:
            self._set_external_forces = self.apply_gravity

        if self._length_scale is None:
            self._length_scale = 1  # 1 => m, 1/100. => cm

        self._out = h5py.File(self._io_filename, self._mode)
        self._data = group(self._out, 'data')
        self._ref = group(self._data, 'ref')
        self._joints = group(self._data, 'joints')
        self._static_data = data(self._data, 'static', 9)
        self._dynamic_data = data(self._data, 'dynamic', 9)
        self._cf_data = data(self._data, 'cf', 15)
        self._solv_data = data(self._data, 'solv', 4)
        self._input = group(self._data, 'input')
        self._nslaws = group(self._data, 'nslaws')

        if self._shape_filename is None:
            self._shape = ShapeCollection(self)
        else:
            self._shape = ShapeCollection(self._shape_filename)

        return self

    def __exit__(self, type_, value, traceback):
        self._out.close()


    def apply_gravity(self, body):
        g = constants.g / self._length_scale
        weight = [0, 0, - body.massValue() * g]
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
        Coordinates and orientations of dynamics objects.
        """
        return self._dynamic_data

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

    def importBRepObject(self, name, position, orientation,
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
            body = body_class(position + orientation, velocity, mass, inertia)
            for contactor in contactors:

                # /!\ shared pointer <-> python ref ...
                shape_instantiated = self._shape.get(contactor.name,
                                                     shape_class, face_class, edge_class)
                self._keep.append(shape_instantiated)
                body.addContactShape(shape_instantiated,
                                     contactor.position,
                                     contactor.orientation,
                                     contactor.group)
            self._set_external_forces(body)

            # add the dynamical system to the non smooth
            # dynamical system
            nsds = self._broadphase.model().nonSmoothDynamicalSystem()
            nsds.insertDynamicalSystem(body)
            nsds.setOSI(body, self._osi)
            nsds.setName(body, str(name))


    def importObject(self, name, position, orientation,
                     velocity, contactors, mass, inertia, body_class, shape_class):

        if body_class is None:
            body_class = BulletDS

        if self._broadphase is not None and 'input' in self._data:
            if mass == 0.:
                # a static object
                rbase = btQuaternion(orientation[1],
                                     orientation[2],
                                     orientation[3],
                                     orientation[0])

                for c in contactors:
                    (w, x, y, z) = c.orientation
                    c_orientation = btQuaternion(x, y, z, w)
                    rc_orientation = mul(rbase, c_orientation)

                    c_origin = btVector3(c.position[0],
                                         c.position[1],
                                         c.position[2])

                    rc_origin = quatRotate(rbase, c_origin)

                    rc_sorigin = btVector3(rc_origin.x() + position[0],
                                           rc_origin.y() + position[1],
                                           rc_origin.z() + position[2])

                    static_cobj = btCollisionObject()
                    static_cobj.setCollisionFlags(
                        btCollisionObject.CF_STATIC_OBJECT)

                    self._static_origins.append(rc_sorigin)

                    self._static_orientations.append(rc_orientation)
                    transform = btTransform(rc_orientation)
                    transform.setOrigin(rc_sorigin)
                    self._static_transforms.append(transform)
                    static_cobj.setWorldTransform(transform)

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
                    bws.setInertia(inertia[0], inertia[1], inertia[2])

                body = body_class(bws,
                                  position + orientation,
                                  velocity,
                                  contactors[0].position,
                                  contactors[0].orientation,
                                  contactors[0].group)

                for contactor in contactors[1:]:
                    shape_id = self._shapeid[contactor.name]

                    body.addCollisionShape(self._shape.get(contactor.name),
                                           contactor.position,
                                           contactor.orientation,
                                           contactor.group)

                # set external forces
                self._set_external_forces(body)

                # add the dynamical system to the non smooth
                # dynamical system
                nsds = self._broadphase.model().nonSmoothDynamicalSystem()
                nsds.insertDynamicalSystem(body)
                nsds.setOSI(body, self._osi)
                nsds.setName(body, str(name))

    def importJoint(self, name):
        if self._broadphase is not None:
            topo = self._broadphase.model().nonSmoothDynamicalSystem().\
                   topology()

            joint_class = getattr(Joints,
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

    def importScene(self, body_class, shape_class, face_class, edge_class):
        """
        Import into the broadphase object all the static and dynamic objects
        from hdf5 file
        """

        for shape_name in self._ref:
            self._shapeid[shape_name] = self._ref[shape_name].attrs['id']
            self._number_of_shapes += 1

        def floatv(v):
            return [float(x) for x in v]

        # import dynamical systems
        if self._broadphase is not None and 'input' in self._data:

            for (name, obj) in self._input.items():
                input_ctrs = [ctr for _n_, ctr in obj.items()]
                mass = obj.attrs['mass']
                position = obj.attrs['position']
                orientation = obj.attrs['orientation']
                velocity = obj.attrs['velocity']
                contactors = [Contactor(ctr.attrs['name'],
                                        int(ctr.attrs['group']),
                                        floatv(ctr.attrs['position']),
                                        floatv(ctr.attrs['orientation']))
                              for ctr in input_ctrs]

                if 'inertia' in obj.attrs:
                    inertia = obj.attrs['inertia']
                else:
                    inertia = None


                if True in ('type' in self.shapes()[ctr.attrs['name']].attrs
                            and 'brep' in self.shapes()[ctr.attrs['name']].attrs['type']
                            for ctr in input_ctrs):
                    # Occ object
                    self.importBRepObject(name, floatv(position), floatv(orientation),
                                          floatv(velocity), contactors, float(mass),
                                          inertia, body_class, shape_class, face_class, edge_class)
                else:
                    # Bullet object
                    self.importObject(name, floatv(position), floatv(orientation),
                                      floatv(velocity), contactors, float(mass),
                                      inertia, body_class, shape_class)

            # import nslaws
            for name in self._nslaws:
                self.importNonSmoothLaw(name)

            for name in self.joints():
                self.importJoint(name)

    def outputStaticObjects(self):
        """
        Outputs positions and orientations of static objects
        """
        time = self._broadphase.model().simulation().nextTime()
        idd = -1
        p = 0
        self._static_data.resize(len(self._static_transforms), 0)

        for transform in self._static_transforms:
            position = transform.getOrigin()
            rotation = transform.getRotation()
            self._static_data[p, :] = \
                [time,
                 idd,
                 position.x(),
                 position.y(),
                 position.z(),
                 rotation.w(),
                 rotation.x(),
                 rotation.y(),
                 rotation.z()]
            idd -= 1
            p += 1

    def outputDynamicObjects(self):
        """
        Outputs positions and orientations of dynamic objects
        """

        current_line = self._dynamic_data.shape[0]

        time = self._broadphase.model().simulation().nextTime()

        positions = self._io.positions(self._broadphase.model())

        self._dynamic_data.resize(current_line + positions.shape[0], 0)

        times = np.empty((positions.shape[0], 1))
        times.fill(time)

        tidd = np.arange(1,
                         positions.shape[0] + 1).reshape(
                             positions.shape[0],
                             1)

        self._dynamic_data[current_line:, :] = np.concatenate((times, tidd,
                                                               positions),
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
                    writer=vtk.vtkXMLPolyDataWriter()
                    writer.SetInputData(reader.GetOutput())
                    writer.SetFileName(tmpf[1])
                    writer.Write()

                    shape_data = str_of_file(tmpf[1])

            else:
                assert os.path.splitext(filename)[-1][1:] == 'vtp'
                shape_data = str_of_file(filename)


            self.addMeshShapeFromString(name, shape_data)

    def addBRepFromString(self, name, shape_data):
        """
        Add a brep contained in a string.
        """
        if name not in self._ref:
            shape = self._ref.create_dataset(name, (1,),
                                             dtype=h5py.new_vlen(str))
            shape[:] = shape_data
            shape.attrs['id'] = self._number_of_shapes
            shape.attrs['type'] = 'brep'
            self._shapeid[name] = shape.attrs['id']
            self._number_of_shapes += 1

    def addContactFromBRep(self, name, brepname, contact_type,
                           index, collision_group=0, associated_shape=None):
        """
        Add contact reference from previously added brep.
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

    def addConvexShape(self, name, points):
        """
        Add a convex shape defined by points.
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

    def addObject(self, name, contactors,
                  position,
                  orientation=[1, 0, 0, 0],
                  velocity=[0, 0, 0, 0, 0, 0],
                  mass=0, inertia=None):
        """
        Add an object.
        Contact detection is defined by a list of contactors.
        The initial position is mandatory : [x, y z].
        If the mass is zero this is a static object.
        """

        if len(orientation) == 2:
            # axis + angle
            axis = orientation[0]
            assert len(axis) == 3
            angle = orientation[1]
            assert type(angle) is float
            n = sin(angle / 2) / np.linalg.norm(axis)

            ori = [cos(angle / 2.), axis[0] * n, axis[1] * n, axis[2] * n]
        else:
            # a given quaternion
            ori = orientation

        if name not in self._input:
            obj = group(self._input, name)
            obj.attrs['mass'] = mass
            obj.attrs['position'] = position
            obj.attrs['orientation'] = ori
            obj.attrs['velocity'] = velocity

            if inertia is not None:
                obj.attrs['inertia'] = inertia

            for num, contactor in enumerate(contactors):
                dat = data(obj, '{0}-{1}'.format(contactor.name, num), 0)
                dat.attrs['name'] = contactor.name
                dat.attrs['group'] = contactor.group
                dat.attrs['position'] = contactor.position
                dat.attrs['orientation'] = contactor.orientation

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
            length_scale=1,
            t0=0,
            T=10,
            h=0.0005,
            multipoints_iterations=True,
            theta=0.50001,
            Newton_max_iter=20,
            set_external_forces=None,
            solver=Numerics.SICONOS_FRICTION_3D_NSGS,
            itermax=100000,
            tolerance=1e-8):
        """
        Run a simulation from inputs in hdf5 file.
        parameters are:
          with_timer : use a timer for log output (default False)
          length_scale : gravity is multiplied by this factor.
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

        """

        from Siconos.Kernel import \
            Model, MoreauJeanOSI, TimeDiscretisation,\
            GenericMechanical, FrictionContact, NewtonImpactFrictionNSL

        from Siconos.Numerics import SICONOS_FRICTION_3D_AlartCurnierNewton

        from Siconos.Mechanics.ContactDetection.Bullet import \
            btConvexHullShape, btCollisionObject, \
            btBoxShape, btQuaternion, btTransform, btConeShape, \
            BulletSpaceFilter, cast_BulletR, \
            BulletWeightedShape, BulletDS, BulletTimeStepping

        if set_external_forces is not None:
            self._set_external_forces = set_external_forces

        if time_stepping is None:
            time_stepping = BulletTimeStepping

        if space_filter is None:
            space_filter = BulletSpaceFilter

        # Model
        #
        model = Model(t0, T)

        # (1) OneStepIntegrators
        joints = list(self.joints())

        self._osi = MoreauJeanOSI(theta)

        # (2) Time discretisation --
        timedisc = TimeDiscretisation(t0, h)

        if len(joints) > 0:
            osnspb = GenericMechanical(SICONOS_FRICTION_3D_AlartCurnierNewton)
        else:
            osnspb = FrictionContact(3, solver)

        osnspb.numericsSolverOptions().iparam[0] = itermax
        osnspb.numericsSolverOptions().dparam[0] = tolerance
        osnspb.setMaxSize(16384)
        osnspb.setMStorageType(1)
        #osnspb.setNumericsVerboseMode(False)

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

        k = 1

        self.importScene(body_class, shape_class, face_class, edge_class)

        model.initialize(simulation)

        self.outputStaticObjects()
        self.outputDynamicObjects()

        while simulation.hasNextEvent():

            print ('step', k)

            log(self._broadphase.buildInteractions, with_timer)\
                (model.currentTime())

            log(simulation.computeOneStep, with_timer)()

            log(self.outputDynamicObjects, with_timer)()

            log(self.outputContactForces, with_timer)()

            log(self.outputSolverInfos, with_timer)()

            log(simulation.nextStep, with_timer)()

            log(self._out.flush)()

            print ('')
            k += 1
