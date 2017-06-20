#!/usr/bin/env python

#
# Multiples contactors with multiples non smooth laws
#


from siconos.mechanics.collision.tools import Contactor
from siconos.io.mechanics_io import Hdf5
from math import pi, sin, cos

# length of first branch
l1 = 60

# length of second branch
l2 = 10

# mass 1
m1 = 1

# mass 2
m2 = 2

# radius of first ball
r1 = 1

# radisu of second ball
r2 = 1

# Creation of the hdf5 file for input/output
with Hdf5() as io:

    #
    io.addPrimitiveShape('Arm1', 'Cylinder', (.3, l1))
    io.addPrimitiveShape('Arm2', 'Cylinder', (.3, l2))
    io.addPrimitiveShape('Mass1', 'Sphere', (r1,))
    io.addPrimitiveShape('Mass2', 'Sphere', (r2,))

    # Definition of the ground shape
    ground_slope = 0.29
    ground_size = (100., 400, 2.)
    io.addPrimitiveShape('Ground', 'Box', ground_size)
    io.addPrimitiveShape('Ground2', 'Box', ground_size)
    block_size = (50., 20., 20.)
    io.addPrimitiveShape('Block', 'Box', block_size)

    # Definition of a non smooth law. As no group ids are specified it
    # is between contactors of group id 0.
    io.addNewtonImpactFrictionNSL('contact', mu=0.3)

    #first branch + first mass the center of gravity is at the center of the
    #Mass1
    ppos = [0, -10., 15]
    io.addObject(
        'arm1', [Contactor('Mass1'),
                 Contactor('Arm1', relative_translation=[0, r1 + l1 / 2., 0])],
        translation=ppos,
        orientation=((1, 0, 0), pi / 2),
        velocity=[0, -5, 0, 0, 0, 0],
        mass=m1)

    io.addJoint('joint2', 'arm1',
                pivot_point=[ppos[0], ppos[1], ppos[2] + l1],
                axis=[1, 0, 0],
                joint_class='PivotJointR')

    # the ground object made with the ground shape. As the mass is
    # not given, it is a static object only involved in contact
    # detection.
    shift = -0.5 * ground_size[2] - sin(ground_slope) * ground_size[1] * 0.5
    io.addObject('ground', [Contactor('Ground')],
                 translation=[0, ground_size[1] * 0.3, shift])
    io.addObject('ground2', [Contactor('Ground2')],
                 translation=[0, 0, 0.],
                 orientation=((1, 0, 0), -ground_slope))
    translat = [0., block_size[2] * 0.5 * sin(ground_slope),
                block_size[2] * 0.5 * cos(ground_slope)]
    io.addObject('block', [Contactor('Block')],
                 translation=translat,
                 orientation=((1, 0, 0), -ground_slope),
                 velocity=[0, 0, 0., 0, 0, 0], mass=m2,)
    # io.addObject('blockref', [Contactor('Block')],
    #              translation=[0., 0., 0.],
    #              orientation=((1, 0, 0), -ground_slope), mass=m2)


# Run the simulation from the inputs previously defined and add
# results to the hdf5 file. The visualisation of the output may be done
# with the vview command.
with Hdf5(mode='r+') as io:

    # By default earth gravity is applied and the units are those
    # of the International System of Units.
    # Because of fixed collision margins used in the collision detection,
    # sizes of small objects may need to be expressed in cm or mm.
    io.run(T=30.)
