#!/usr/bin/env python

#
# Multiples contactors with multiples non smooth laws
#

from siconos.mechanics.collision.tools import Contactor
from siconos.io.mechanics_io import Hdf5
import siconos.io.mechanics_io

# Creation of the hdf5 file for input/output
with Hdf5(mode='w') as io:

    # Definition of two boxes shape
    io.add_primitive_shape('BigBox', 'Box', (3, 5, 2))

    io.add_primitive_shape('LongBox', 'Box', (1, 1, 5))

    # Definition of the ground shape
    io.add_primitive_shape('Ground', 'Box', (100, 100, .5), insideMargin=0.04)

    # Definition of a non smooth law between groups 0 and 1
    io.add_Newton_impact_friction_nsl('contact1', mu=0.1,
                                  collision_group1=0, collision_group2=1)

    # Definition of a non smooth law between groups 1 and 1
    io.add_Newton_impact_friction_nsl('contact2', mu=0.7,
                                  collision_group1=1, collision_group2=1)

    # Definition of a non smooth law between groups 0 and 0
    io.add_Newton_impact_friction_nsl('contact3', mu=0.1,
                                  collision_group1=0, collision_group2=0)

    # A 'two boxes object made with two Contactors.
    # As a mass is given, it is a dynamic system involved in contact
    # detection and in the simulation.
    io.add_object('twoboxes', [Contactor('BigBox',
                                        collision_group=0,
                                        relative_translation=[0, 0, 0]),
                                 Contactor('LongBox', collision_group=1,
                                           relative_translation=[0, 0, 0])],
                 translation=[0, 0, 3],
                 velocity=[10, 0, 0, 1, 1, 1],
                 mass=1)

    # the ground object made with the ground shape. As the mass is
    # not given, it is a static object only involved for contact
    # detection.
    io.add_object('ground', [Contactor('Ground', collision_group=1)],
                 translation=[0, 0, 0])


# Run the simulation from the inputs previously defined and add
# results to the hdf5 file. The visualisation of the output may be done
# with the vview command.
with Hdf5(mode='r+') as io:

    # By default earth gravity is applied and the units are those
    # of the International System of Units.
    # Because of fixed collision margins used in the collision detection,
    # sizes of small objects may need to be expressed in cm or mm.
    io.run()
