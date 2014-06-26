# simple bullet input and output

from __future__ import print_function

import os
import sys

from math import cos, sin

import VtkShapes
import shlex
import numpy as np
import h5py

from Siconos.Mechanics.ContactDetection import Contactor

from Siconos.Mechanics.ContactDetection.Bullet import \
    BulletDS, BulletWeightedShape, \
    btCollisionObject, btQuaternion, btTransform, btVector3, quatRotate

from Siconos.Mechanics.ContactDetection.Bullet import \
    cast_BulletR

from Siconos.Mechanics.ContactDetection.Bullet.BulletWrap import __mul__ as mul


from Siconos.Kernel import cast_NewtonImpactFrictionNSL, EqualityConditionNSL, \
    Interaction

import Siconos.Kernel as Kernel
import Siconos.Mechanics as Mechanics

from Siconos.IO import MechanicsIO

import Siconos.Numerics as Numerics

from scipy import constants

import time

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
        return f.read()

class Dat():
    """a Dat context manager reads at instantiation the positions and
       orientations of collision objects from :

       - a ref file (default ref.txt) with shape primitives or shape
         url

       - an input .dat file (default is input.dat)

       input format is :
       shaped_id object_group mass px py pz ow ox oy oz vx vy vx vo1 vo2 vo3

       with
         shape_id : line number in ref file (an integer)
         object group : an integer ; negative means a static object
         mass : mass of the object (a float)
         px py pz : position (float)
         ow ox oy oz : orientation (as an unit quaternion)
         vx vy vx vo1 vo2 vo3 : velocity

       It provides functions to output position and orientation during
       simulation (output is done by default in pos.dat)

       output format is : time object_id px py pz ow ox oy oz

       with:
         time : float
         object_id : the object id (int)
         px, py, pz : components of the position (float)
         ow, ox, oy oz : components of an unit quaternion (float)
    """

    def __init__(self, broadphase, osi, shape_filename='ref.txt',
                 input_filename='input.dat',
                 set_external_forces=apply_gravity):
        self._broadphase = broadphase
        self._osi = osi
        self._input_filename = input_filename
        self._static_origins = []
        self._static_orientations = []
        self._static_transforms = []
        self._static_cobjs = []
        self._shape = VtkShapes.Collection(shape_filename)
        self._static_pos_file = None
        self._dynamic_pos_file = None
        self._contact_forces_file = None
        self._solver_traces_file = None
        self._io = MechanicsIO()

        # read data
        with open(self._input_filename, 'r') as input_file:

            with open('bindings.dat', 'w') as bind_file:

                ids = -1
                idd = 1
                for line in input_file:
                    sline = shlex.split(line)
                    if len(sline) > 3:
                        shape_id = int(sline[0])
                        group_id = int(sline[1])
                        mass = float(sline[2])
                        q0, q1, q2, w, x, y, z, v0, v1, v2, v3, v4, v5 =\
                          [ float(i) for i in sline[3:]]

                        if group_id < 0:
                            # a static object
                            static_cobj = btCollisionObject()
                            static_cobj.setCollisionFlags(
                                btCollisionObject.CF_STATIC_OBJECT)
                            origin = btVector3(q0, q1, q2)
                            self._static_origins.append(origin)
                            orientation = btQuaternion(x, y, z, w)
                            self._static_orientations.append(orientation)
                            transform = btTransform(orientation)
                            transform.setOrigin(origin)
                            self._static_transforms.append(transform)
                            static_cobj.setWorldTransform(transform)
                            static_cobj.setCollisionShape(
                                self._shape.at_index(shape_id))
                            self._static_cobjs.append(static_cobj)
                            broadphase.addStaticObject(static_cobj, abs(group_id)-1)
                            bind_file.write('{0} {1}\n'.format(ids, shape_id))
                            ids -= 1

                        else:
                            # a moving object
                            body = BulletDS(BulletWeightedShape(
                                self._shape.at_index(shape_id), mass),
                            [q0, q1, q2, w, x, y, z],
                            [v0, v1, v2, v3, v4, v5])

                              # set external forces
                            set_external_forces(body)

                            # add the dynamical system to the non smooth
                            # dynamical system
                            self._broadphase.model().nonSmoothDynamicalSystem().\
                            insertDynamicalSystem(body)
                            self._osi.insertDynamicalSystem(body)
                            bind_file.write('{0} {1}\n'.format(idd, shape_id))
                            idd += 1

    def __enter__(self):
        self._static_pos_file = open('spos.dat', 'w')
        self._dynamic_pos_file = open('dpos.dat', 'w')
        self._contact_forces_file = open('cf.dat', 'w')
        self._solver_traces_file = open('solv.dat', 'w')
        return self

    def __exit__(self, type_, value, traceback):
        self._contact_forces_file.close()
        self._static_pos_file.close()
        self._dynamic_pos_file.close()
        self._solver_traces_file.close()


    def outputStaticObjects(self):
        """
        Outputs positions and orientations of static objects
        """
        time = self._broadphase.model().simulation().nextTime()
        idd = -1
        for transform in self._static_transforms:
            position = transform.getOrigin()
            rotation = transform.getRotation()
            self._static_pos_file.write(
                '{0} {1} {2} {3} {4} {5} {6} {7} {8}\n'.
                format(time,
                       idd,
                       position.x(),
                       position.y(),
                       position.z(),
                       rotation.w(),
                       rotation.x(),
                       rotation.y(),
                       rotation.z()))
            idd -= 1

    def outputDynamicObjects(self):
        """
        Outputs positions and orientations of dynamic objects
        """
        time = self._broadphase.model().simulation().nextTime()

        positions = self._io.positions(self._broadphase.model())

        times = np.empty((positions.shape[0], 1))
        times.fill(time)

        tidd = np.arange(1,
                         positions.shape[0] + 1).reshape(positions.shape[0], 1)

        np.savetxt(self._dynamic_pos_file,
                   np.concatenate((times, tidd, positions),
                                  axis=1))

    def outputContactForces(self):
        """
        Outputs contact forces
        """
        if self._broadphase.model().nonSmoothDynamicalSystem().\
                topology().indexSetsSize() > 1:
            time = self._broadphase.model().simulation().nextTime()
            for inter in self._broadphase.model().nonSmoothDynamicalSystem().\
                topology().indexSet(1).interactions():
                bullet_relation = cast_BulletR(inter.relation())
                if bullet_relation is not None:
                    nslaw = inter.nslaw()
                    mu = cast_NewtonImpactFrictionNSL(nslaw).mu()
                    nc = bullet_relation.nc()
                    lambda_ = inter.lambda_(1)
                    if(True):
                        jachqt = bullet_relation.jachqT()
                        cf = np.dot(jachqt.transpose(), lambda_)
                        cp = bullet_relation.contactPoint()
                        posa = cp.getPositionWorldOnA()
                        posb = cp.getPositionWorldOnB()
                        self._contact_forces_file.write(
                        '{0} {1} {2} {3} {4} {5} {6} {7} {8} {9} {10} {11} {12} {13}\n'.
                        format(time,
                               mu,
                               posa.x(),
                               posa.y(),
                               posa.z(),
                               posb.x(),
                               posb.y(),
                               posb.z(),
                               nc[0], nc[1], nc[2],
                               cf[0], cf[1], cf[2]))

    def outputSolverInfos(self):
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

        self._solver_traces_file.write('{0} {1} {2} {3}\n'.
                                       format(time, iterations, precision,
                                              local_precision))


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
                 set_external_forces=apply_gravity):

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

    def __enter__(self):
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
            self._shape = VtkShapes.Collection(self._out)
        else:
            self._shape = VtkShapes.Collection(self._shape_filename)

        return self

    def __exit__(self, type_, value, traceback):
        self._out.close()


    def shapes(self):
        return self._ref

    def static_data(self):
        return self._static_data

    def dynamic_data(self):
        return self._dynamic_data

    def contact_forces_data(self):
        return self._cf_data

    def solver_data(self):
        return self._solv_data

    def instances(self):
        return self._input

    def nonsmooth_laws(self):
        return self._nslaws

    def joints(self):
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

    def importObject(self, name, position, orientation,
                     velocity, contactors, mass):

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
                    shape_id = self._shapeid[c.name]
                    static_cobj.setCollisionShape(
                        self._shape.at_index(shape_id))
                    self._static_cobjs.append(static_cobj)
                    self._broadphase.addStaticObject(static_cobj,
                                                     int(c.group))

            else:
                # a moving object
                shape_id = self._shapeid[contactors[0].name]
                body = BulletDS(BulletWeightedShape(
                    self._shape.at_index(shape_id), mass),
                    position + orientation,
                    velocity,
                    contactors[0].position,
                    contactors[0].orientation,
                    contactors[0].group)

                for contactor in contactors[1:]:
                    shape_id = self._shapeid[contactor.name]

                    body.addCollisionShape(self._shape.at_index(shape_id),
                                           contactor.position,
                                           contactor.orientation,
                                           contactor.group)

                # set external forces
                self._set_external_forces(body)

                # add the dynamical system to the non smooth
                # dynamical system
                self._broadphase.model().nonSmoothDynamicalSystem().\
                     insertDynamicalSystem(body, str(name))
                self._osi.insertDynamicalSystem(body)

    def importJoint(self, name):
        if self._broadphase is not None:
            topo = self._broadphase.model().nonSmoothDynamicalSystem().\
                   topology()
            joint_class = getattr(Mechanics.Joints,
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

    def importScene(self):
        """
        Import into the broadphase object all the static and dynamic objects
        from hdf5 file
        """

        for shape_name in self._ref:
            self._shapeid[shape_name] = self._ref[shape_name].attrs['id']
            self._number_of_shapes += 1

        self._shape = VtkShapes.Collection(self._out)

        def floatv(v):
            return [float(x) for x in v]

        # import dynamical systems
        if self._broadphase is not None and 'input' in self._data:

            for (name, obj) in self._input.items():
                mass = obj.attrs['mass']
                position = obj.attrs['position']
                orientation = obj.attrs['orientation']
                velocity = obj.attrs['velocity']
                contactors = [Contactor(ctr.attrs['name'],
                                        int(ctr.attrs['group']),
                                        floatv(ctr.attrs['position']),
                                        floatv(ctr.attrs['orientation']))
                              for _n_, ctr in obj.items()]

                self.importObject(name, floatv(position), floatv(orientation),
                                  floatv(velocity), contactors, float(mass))

            # import nslaws
            for name in self._nslaws:
                self.importNonSmoothLaw(name)

            for name in self.joints():
                print ('import joint:', name)
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

    def insertShapeFromFile(self, name, filename):
        """
        Insert a mesh shape from a file.
        Accepted format : mesh encoded in VTK .vtp file
        """
        if name not in self._ref:
            shape = self._ref.create_dataset(name, (1,),
                                             dtype=h5py.new_vlen(str))
            shape[:] = str_of_file(filename)
            shape.attrs['id'] = self._number_of_shapes
            self._shapeid[name] = shape.attrs['id']
            self._shape = VtkShapes.Collection(self._out)
            self._number_of_shapes += 1

    def insertConvexShape(self, name, points):
        """
        Insert a convex shape defined by points.
        """
        if name not in self._ref:
            apoints = np.array(points)
            shape = self._ref.create_dataset(name,
                                             (apoints.shape[0],
                                              apoints.shape[1]))
            shape[:] = points[:]
            shape.attrs['id'] = self._number_of_shapes
            self._shapeid[name] = shape.attrs['id']
            self._shape = VtkShapes.Collection(self._out)
            self._number_of_shapes += 1

    def insertPrimitiveShape(self, name, primitive, params):
        """
        Insert a primitive shape.
        """
        if name not in self._ref:
            shape = self._ref.create_dataset(name, (1, len(params)))
            shape.attrs['id'] = self._number_of_shapes
            shape.attrs['primitive'] = primitive
            shape[:] = params
            self._shapeid[name] = shape.attrs['id']
            self._shape = VtkShapes.Collection(self._out)
            self._number_of_shapes += 1

    def insertObject(self, name, contactors,
                     position,
                     orientation=[1, 0, 0, 0],
                     velocity=[0, 0, 0, 0, 0, 0],
                     mass=0):
        """
        Insertion of an object.
        Contact detection is defined by a list of contactors.
        The initial position is mandatory : [x, y z]
        if the mass is zero this a static object.
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

    def insertNewtonImpactFrictionNSL(self, name, mu, e=0, collision_group1=0,
                                      collision_group2=0):
        """
        Insert a nonsmooth law for contact between 2 groups.
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


    def insertJoint(self, name, object1, object2=None, pivot_point=[0, 0, 0],
                    axis=[0, 1, 0],
                    joint_class='PivotJointR'):
        """
        insert a pivot joint between two objects
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
            t0=0,
            T=10,
            h=0.0005,
            multipointIterations=True,
            theta=0.50001,
            set_external_forces=None,
            solver=Numerics.SICONOS_FRICTION_3D_NSGS,
            NewtonMaxIter=20,
            itermax=100000,
            tolerance=1e-8):

        from Siconos.Kernel import \
            Model, MoreauJeanOSI, TimeDiscretisation,\
            GenericMechanical, FrictionContact, NewtonImpactFrictionNSL,\
            BlockCSRMatrix

        from Siconos.Numerics import SICONOS_FRICTION_3D_AlartCurnierNewton

        from Siconos.Mechanics.ContactDetection.Bullet import IO, \
            btConvexHullShape, btCollisionObject, \
            btBoxShape, btQuaternion, btTransform, btConeShape, \
            BulletSpaceFilter, cast_BulletR, \
            BulletWeightedShape, BulletDS, BulletTimeStepping

        if set_external_forces is not None:
            self._set_external_forces = set_external_forces

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
        self._broadphase = BulletSpaceFilter(model)
        if not multipointIterations:
            print("""
            ConvexConvexMultipointIterations and PlaneConvexMultipointIterations are unset
            """)
        else:
            self._broadphase.collisionConfiguration().\
                setConvexConvexMultipointIterations()
            self._broadphase.collisionConfiguration().\
                setPlaneConvexMultipointIterations()

        # (6) Simulation setup with (1) (2) (3) (4) (5)
        simulation = BulletTimeStepping(timedisc, self._broadphase)
        simulation.insertIntegrator(self._osi)
        simulation.insertNonSmoothProblem(osnspb)
        simulation.setNewtonMaxIteration(NewtonMaxIter)

        k = 1

        self.importScene()

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
