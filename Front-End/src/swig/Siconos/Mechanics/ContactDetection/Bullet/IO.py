# simple bullet input and output

import VtkShapes
import shlex
import numpy as np

from Siconos.Mechanics.ContactDetection.Bullet import \
    BulletDS, BulletWeightedShape, \
    btCollisionObject, btQuaternion, btTransform, btVector3

from Siconos.Mechanics.ContactDetection.Bullet import \
    cast_BulletR

from Siconos.Kernel import cast_NewtonImpactFrictionNSL

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
                            self._reference[ids] = shape_id
                            self._static_cobjs.append(static_cobj)
                            broadphase.addStaticObject(static_cobj)
                            broadphase.addStaticShape(self._shape.at_index(shape_id))
                            bind_file.write('{0} {1}\n'.format(ids, shape_id))
                            ids -= 1

                        else:
                            # a moving object
                            body = BulletDS(BulletWeightedShape(
                                self._shape.at_index(shape_id), mass),
                            [q0, q1, q2, w, x, y, z],
                            [v0, v1, v2, v3, v4, v5])
                            self._reference[idd] = \
                              shape_id

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
                topology().indexSet(1).vertices():
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
        so = self._broadphase.model().simulation().oneStepNSProblem(0).numericsSolverOptions()
        iterations = so.iparam[3]
        precision = so.dparam[2]
        local_precision = so.dparam[3]
        self._solver_traces_file.write('{0} {1} {2} {3}\n'.format(time, iterations, precision, local_precision))
