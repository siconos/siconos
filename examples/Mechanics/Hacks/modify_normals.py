#!/usr/bin/env python

#
# Example of overriding BulletR to modify the contact normals.
#

from siconos.mechanics.collision.tools import Contactor
from siconos.io.mechanics_io import Hdf5
import siconos.numerics as Numerics
import siconos.mechanics as Mechanics
import numpy as np

# Creation of the hdf5 file for input/output
with Hdf5() as io:

    # Definition of a cube as a convex shape
    io.addPrimitiveShape('Cube', 'Box', (1,1,1))

    # Definition of the ground shape
    io.addPrimitiveShape('Ground', 'Box', (100, 100, .5))

    # Definition of a non smooth law. As no group ids are specified it
    # is between contactors of group id 0.
    io.addNewtonImpactFrictionNSL('contact', mu=0.3)

    # The cube object made with an unique Contactor : the cube shape.
    # As a mass is given, it is a dynamic system involved in contact
    # detection and in the simulation.  With no group id specified the
    # Contactor belongs to group 0
    io.addObject('cube', [Contactor('Cube')], translation=[0, 0, 2],
                 velocity=[10, 0, 0, 1, 1, 1],
                 mass=1)

    # the ground object made with the ground shape. As the mass is
    # not given, it is a static object only involved in contact
    # detection.
    io.addObject('ground', [Contactor('Ground')],
                 translation=[0, 0, 0])

# Need to keep the MyBulletR relations somewhere otherwise some
# SWIG-related refcounting problems cause a crash.
relations_keeper = {}

class MyBulletR(Mechanics.collision.bullet.BulletR):
    def updateContactPoints(self, point):
        # Call usual updateContactPoints
        super(self.__class__,self).updateContactPoints(point)

        # Add some noise to the normal's direction
        n = self.nc() + np.random.normal(0, 0.1, 3)
        n = n / np.linalg.norm(n)
        self.setnc(n)

    # This function is called when the contact Interaction is no
    # longer needed.  We can remove it from the keeper.
    # TODO: There appears to be a left-over refcount somewhere, memory
    # is not correctly freed here and will accumulate.
    def preDelete(self):
        del relations_keeper[self]

class MyBulletManager(Mechanics.collision.bullet.SiconosBulletCollisionManager):
    # Call the appropriate constructor depending on whether options
    # were provided -- it appears that SWIG is not smart enough to do
    # this, and will happily try to dereference a null pointer
    # otherwise.
    def __init__(self, model, options):
        if options: super(self.__class__,self).__init__(options)
        else:  super(self.__class__,self).__init__()

    # Override makeBulletR to return an instance of our own
    # BulletR-derived class.
    def makeBulletR(self, ds1, shape1, ds2, shape2, manifoldpoint,
                    flip=False, y_correction_A=0, y_correction_B=0, scaling=1):
        r = MyBulletR(manifoldpoint, flip, y_correction_A, y_correction_B, scaling)
        relations_keeper[r] = True
        return r

# Run the simulation from the inputs previously defined and add
# results to the hdf5 file. The visualisation of the output may be done
# with the vview command.
with Hdf5(mode='r+') as io:

    # We assign our custom CollisionManager to the "space_filter" input
    io.run(with_timer=False,
            time_stepping=None,
            space_filter=MyBulletManager,
            body_class=None,
            shape_class=None,
            face_class=None,
            edge_class=None,
            gravity_scale=1,
            t0=0,
            T=10,
            h=0.005,
            theta=0.50001,
            Newton_max_iter=20,
            set_external_forces=None,
            solver=Numerics.SICONOS_FRICTION_3D_NSGS,
            itermax=100000,
            tolerance=1e-8,
            numerics_verbose=False,
            output_frequency=None)
