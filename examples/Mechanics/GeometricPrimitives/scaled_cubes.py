#!/usr/bin/env python

#
# Here we demonstrate the use of two very small cubes simulated with
# geometry at 1000x its scale.  Without scaling, contact detection is
# incorrect for an object of this size.  When contact detection is
# performed with scaled-up geometries, the contact points are
# correctly generated and used by Siconos.
#

from siconos.mechanics.collision.tools import Contactor
from siconos.io.mechanics_run import MechanicsHdf5Runner
import siconos.numerics as Numerics

# We need to pass some options to the Bullet backend
from siconos.mechanics.collision.bullet import SiconosBulletOptions
options = SiconosBulletOptions()
options.worldScale = 1000.0

# Creation of the hdf5 file for input/output
with MechanicsHdf5Runner() as io:

    # Definition of a cube as a convex shape
    io.add_convex_shape('CubeCH', [
        (-0.001, 0.001, -0.001),
        (-0.001, -0.001, -0.001),
        (-0.001, -0.001, 0.001),
        (-0.001, 0.001, 0.001),
        (0.001, 0.001, 0.001),
        (0.001, 0.001, -0.001),
        (0.001, -0.001, -0.001),
        (0.001, -0.001, 0.001)])

    # Definition of a cube as a primitive shape
    io.add_primitive_shape('CubeP', 'Box', (0.002, 0.002, 0.002))

    # Definition of the ground shape
    io.add_primitive_shape('Ground', 'Box', (0.1, 0.1, .01))

    # Definition of a non smooth law. As no group ids are specified it
    # is between contactors of group id 0.
    io.add_Newton_impact_friction_nsl('contact', mu=0.1)

    # The cube object made with an unique Contactor : the cube shape.
    # As a mass is given, it is a dynamic system involved in contact
    # detection and in the simulation.  With no group id specified the
    # Contactor belongs to group 0
    io.add_object('cubeCH', [Contactor('CubeCH')],
                  translation=[0, 0.003, 0.005],
                  velocity=[0.1, 0, 0, 1, 1, 1],
                  mass=0.1)

    # The primitive cube geometry object
    io.add_object('cubeP', [Contactor('CubeP')],
                  translation=[0, -0.003, 0.005],
                  velocity=[0.1, 0, 0, 1, 1, 1],
                  mass=0.1)

    # the ground object made with the ground shape. As the mass is
    # not given, it is a static object only involved in contact
    # detection.
    io.add_object('ground', [Contactor('Ground')],
                  translation=[0, 0, -0.005])


# Run the simulation from the inputs previously defined and add
# results to the hdf5 file. The visualisation of the output may be done
# with the vview command.
with MechanicsHdf5Runner(mode='r+') as io:

    # By default earth gravity is applied and the units are those
    # of the International System of Units.
    # Because of fixed collision margins used in the collision detection,
    # sizes of small objects may need to be expressed in cm or mm.
    io.run(with_timer=False,
           options=options,
           t0=0,
           T=10,
           h=0.005,
           theta=0.50001,
           Newton_max_iter=4,
           set_external_forces=None,
           solver=Numerics.SICONOS_FRICTION_3D_NSGS,
           itermax=10000,
           tolerance=1e-8,
           numerics_verbose=False,
           output_frequency=None)
