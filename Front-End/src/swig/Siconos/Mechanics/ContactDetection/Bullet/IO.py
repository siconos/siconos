# simple bullet input and output

import VtkShapes
import shlex
from itertools import chain
import numpy as np

from Siconos.Mechanics.ContactDetection.Bullet import \
    BulletDS, BulletWeightedShape, \
    btCollisionObject, btQuaternion, btTransform

from Siconos.Mechanics.ContactDetection.Bullet import \
    cast_BulletDS, cast_BulletR

from Siconos.IO import MechanicsIO


def object_id(obj):
    """returns an unique object identifier"""
    return obj.__hash__()


def apply_gravity(body):
    g = 9.81
    weight = [0, 0, - body.massValue() * g]
    body.setFExtPtr(weight)


class Dat():
    """a Dat object reads at instantiation the positions and
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
        self._reference = dict()
        self._static_cobjs = []
        self._shape = VtkShapes.Collection(shape_filename)
        self._pos_file = None
        self._contact_forces_file = None
        self._io = MechanicsIO()

        # read data
        with open(self._input_filename, 'r') as input_file:

            for line in input_file:
                sline = shlex.split(line)
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
                    static_cobj.setWorldTransform(btTransform(
                        btQuaternion(x, y, z, w)))
                    static_cobj.setCollisionShape(
                        self._shape.at_index(shape_id))
                    self._reference[object_id(static_cobj)] = shape_id
                    self._static_cobjs.append(static_cobj)
                    broadphase.addStaticObject(static_cobj)
                    broadphase.addStaticShape(self._shape.at_index(shape_id))
                else:
                    # a moving object
                    body = BulletDS(BulletWeightedShape(
                        self._shape.at_index(shape_id), mass),
                        [q0, q1, q2, w, x, y, z],
                        [v0, v1, v2, v3, v4, v5])
                    self._reference[object_id(body.collisionObject())] = \
                        shape_id

                    # set external forces
                    set_external_forces(body)

                    # add the dynamical system to the non smooth
                    # dynamical system
                    self._broadphase.model().nonSmoothDynamicalSystem().\
                        insertDynamicalSystem(body)
                    self._osi.insertDynamicalSystem(body)

    def __enter__(self):
        self._pos_file = open('pos.dat', 'w')
        self._contact_forces_file = open('contact_forces.dat', 'w')
        return self

    def __exit__(self, type_, value, traceback):
        with open('bindings.dat', 'w') as bind_file:
            for collision_object in chain(
                self._broadphase.staticObjects(),
                    [cast_BulletDS(ds).collisionObject()
                     for ds in self._broadphase.model().
                     nonSmoothDynamicalSystem().topology().dSG(0).vertices()]):
                bind_file.write('{0} {1}\n'.format(object_id(collision_object),
                                self._reference[object_id(collision_object)]))

        self._contact_forces_file.close()
        self._pos_file.close()

    def outputStaticObjects(self):
        """
        Outputs positions and orientations of static objects
        """
        time = self._broadphase.model().simulation().nextTime()
        for collision_object in self._broadphase.staticObjects():
            position = collision_object.getWorldTransform().getOrigin()
            rotation = collision_object.getWorldTransform().getRotation()
            self._pos_file.write('{0} {1} {2} {3} {4} {5} {6} {7} {8}\n'.
                                 format(time, object_id(collision_object),
                                        position.x(),
                                        position.y(),
                                        position.z(),
                                        rotation.w(),
                                        rotation.x(),
                                        rotation.y(),
                                        rotation.z()))

    def outputDynamicObjects(self):
        """
        Outputs positions and orientations of dynamic objects
        """
        time = self._broadphase.model().simulation().nextTime()

        ids = self._io.dynamicIds(self._broadphase.model())
        positions = self._io.positions(self._broadphase.model())

        id_= 0
        for row in positions:
            obj_id = ids[id_]
            id_ += 1
            self._pos_file.write('{0} {1} {2} {3} {4} {5} {6} {7} {8}\n'.
                                 format(time, obj_id,
                                        row[0], row[1], row[2], 
                                        row[3], row[4], row[5], row[6]))

    def outputContactForces(self):
        """
        Outputs contact forces
        """
        if self._broadphase.model().nonSmoothDynamicalSystem().\
                topology().indexSetsSize() > 1:
            time = self._broadphase.model().simulation().nextTime()
            for inter in self._broadphase.model().nonSmoothDynamicalSystem().\
                  topology().indexSet(1).vertices():
                bullet_relation = cast_BulletR(inter.relation())
                # nslaw = inter.nslaw()
                # cast nslaw if NewtonImpactFrictionNSL [...]
                nc = bullet_relation.nc()
                lambda_ = inter.lambda_(1)
                if not (lambda_[0] == 0. and lambda_[1] == 0. 
                        and lambda_[2] == 0.):
                    jachqt = bullet_relation.jachqT()
                    cf = np.dot(jachqt.transpose(), lambda_)
                    cm = bullet_relation.contactManifold()
                    if (cm.getNumContacts() > bullet_relation.contactNum()):
                        cp = cm.getContactPoint(bullet_relation.contactNum())
                        posa = cp.getPositionWorldOnA()
                        self._contact_forces_file.write(
                        '{0} {1} {2} {3} {4} {5} {6} {7} {8} {9}\n'.
                        format(time,
                               posa.x(),
                               posa.y(),
                               posa.z(),
                               nc[0], nc[1], nc[2],
                               cf[0], cf[1], cf[2]))
