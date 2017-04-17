#!/usr/bin/env python

#
# Example of a double pendulum
#

from siconos.mechanics.collision.tools import Contactor
from siconos.io.mechanics_io import Hdf5

from math import pi

# length of first branch
l1 = 30

# length of second branch
l2 = 10

# mass 1
m1 = 1

# mass 2
m2 = 10

# radius of first ball
r1 = 1

# radisu of second ball
r2 = 1

# gap between ball branch to avoid friction
gap=0.

# gap between ground
hgap = 0.1

# size of a brick
bx = 5
by = 2
bz = 2

# Creation of the hdf5 file for input/output
with Hdf5() as io:

    #
    io.addPrimitiveShape('Arm1', 'Cylinder', (.3, l1))
    io.addPrimitiveShape('Arm2', 'Cylinder', (.3, l2))
    io.addPrimitiveShape('Mass1', 'Sphere', (r1,))
    io.addPrimitiveShape('Mass2', 'Sphere', (r2,))

    # Definition of the ground shape
    io.addPrimitiveShape('Ground', 'Box', (100, 100, .5))

    # the brick shape
    io.addPrimitiveShape('Brick', 'Box', (bx, by, bz))

    # Definition of a non smooth law. As no group ids are specified it
    # is between contactors of group id 0.
    io.addNewtonImpactFrictionNSL('contact', mu=0.1)

    # first branch + first mass the center of gravity is at the center of the
    # Mass1
    io.addObject('arm1', [Contactor('Mass1'),
                          Contactor('Arm1',
                                    relative_translation=[0, r1+l1/2., 0])],
                    translation=[0, 0, r2 + gap + r2 + l2 + r1 + hgap],
                    orientation=((1, 0, 0), pi/2),
                    mass=m1)

    # second branch + second mass
    io.addObject('arm2', [Contactor('Mass2'),
                          Contactor('Arm2',
                                    relative_translation=[0, r2+l2/2., 0])],
                    translation=[0, 0, r2 + gap],
                    orientation=((1, 0, 0), pi/2),
                    velocity=[0, 20, 0, 0, 0, 0],
                    mass=m2)

    io.addJoint('joint1', 'arm1', 'arm2', [0, -r1, 0],
                [1, 0, 0],
                'PivotJointR')

    io.addJoint('joint2', 'arm1',
                pivot_point=[0, 0, r2 + gap + r2 + l2 + r1 + hgap + l1],
                axis=[1, 0, 0],
                joint_class='PivotJointR')

    # a brick wall
    H = 3   # heigh
    L = 2   # length
    for k in range(0, H-1):
        for n in range(0, L):
            io.addObject('brick{0}'.format(k+n*H),
                            [Contactor('Brick')],
                            translation=[n*bx-L*bx/2. + (k % 2) * bx/2.,
                                         -5,
                                         k*bz + bz/2.], mass=2)

    k = H-1
    for n in range(1, L):
        io.addObject('brick{0}'.format(k+n*H),
                        [Contactor('Brick')],
                        translation=[n*bx-L*bx/2. + (k % 2) * bx/2.,
                                     -5,
                                     k*bz + bz/2.], mass=2)

    # the ground object made with the ground shape. As the mass is
    # not given, it is a static object only involved in contact
    # detection.
    io.addObject('ground', [Contactor('Ground')],
                 translation=[0, 0, -.25])


# Run the simulation from the inputs previously defined and add
# results to the hdf5 file. The visualisation of the output may be done
# with the vview command.
with Hdf5(mode='r+') as io:

    io.run()
