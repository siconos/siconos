#!/usr/bin/env python

#
# Example of one object under gravity with one contactor and a ground
#

from siconos.mechanics.collision.tools import Contactor
from siconos.io.mechanics_io import Hdf5
import siconos.numerics as Numerics
import siconos.kernel as Kernel

# Creation of the hdf5 file for input/output
with Hdf5() as io:

    # Definition of a cube as a convex shape
    io.addConvexShape('Cube', [
        (-1.0, 1.0, -1.0),
        (-1.0, -1.0, -1.0),
        (-1.0, -1.0, 1.0),
        (-1.0, 1.0, 1.0),
        (1.0, 1.0, 1.0),
        (1.0, 1.0, -1.0),
        (1.0, -1.0, -1.0),
        (1.0, -1.0, 1.0)])

    # Alternative to the previous convex shape definition.
    # io.addPrimitiveShape('Cube1', 'Box', (2, 2, 2))

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


# Run the simulation from the inputs previously defined and add
# results to the hdf5 file. The visualisation of the output may be done
# with the vview command.
with Hdf5(mode='r+') as io:

    # By default earth gravity is applied and the units are those
    # of the International System of Units.
    # Because of fixed collision margins used in the collision detection,
    # sizes of small objects may need to be expressed in cm or mm.
    io.run(with_timer=False,
            time_stepping=Kernel.TimeSteppingDirectProjection,
            osi=Kernel.MoreauJeanDirectProjectionOSI,
            space_filter=None,
            body_class=None,
            shape_class=None,
            face_class=None,
            edge_class=None,
            gravity_scale=1,
            t0=0,
            T=10,
            h=0.1,
            theta=0.50001,
            Newton_max_iter=20,
            set_external_forces=None,
            solver=Numerics.SICONOS_FRICTION_3D_NSGS,
            itermax=100000,
            tolerance=1e-8,
            numerics_verbose=False,
            output_frequency=None,
            Newton_update_interactions=True,
            projection_itermax=5,
            projection_tolerance=1e-8,
            projection_tolerance_unilateral=1e-8,
    )
