#!/usr/bin/env python

#
# Example of a double pendulum
#

from siconos.mechanics.collision.tools import Contactor
from siconos.io.mechanics_run import MechanicsHdf5Runner
import siconos
from OCC.BRepPrimAPI import BRepPrimAPI_MakeBox, BRepPrimAPI_MakeSphere, BRepPrimAPI_MakeCylinder
from OCC.gp import gp_Pnt

from math import pi
import numpy as np

siconos.io.mechanics_run.set_backend('occ')


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

# gap
gap = .47

# hgap
hgap = 0.

# size of a brick
bx = 5
by = 2
bz = 2

sphere1 = BRepPrimAPI_MakeSphere(r1).Shape()
sphere2 = BRepPrimAPI_MakeSphere(r2).Shape()
cylinder1 = BRepPrimAPI_MakeCylinder(.3, l1).Shape()
cylinder2 = BRepPrimAPI_MakeCylinder(.3, l2).Shape()
ground = BRepPrimAPI_MakeBox(gp_Pnt(-50, -50, 0), 100., 100., .5).Shape()
brick = BRepPrimAPI_MakeBox(gp_Pnt(-bx/2.,-by/2.,-bz/2.),bx, by, bz).Shape()

# Creation of the hdf5 file for input/output
with MechanicsHdf5Runner() as io:

    #
    io.add_occ_shape('Mass1', sphere1)
    io.add_occ_shape('Mass2', sphere2)
    io.add_occ_shape('Arm1', cylinder1)
    io.add_occ_shape('Arm2', cylinder2)
    io.add_occ_shape('Ground', ground)
    io.add_occ_shape('Brick', sphere1)

    # Definition of a non smooth law. As no group ids are specified it
    # is between contactors of group id 0.
    io.add_Newton_impact_friction_nsl('contact', mu=0.1)

    # first branch + first mass the center of gravity is at the center of the
    # Mass1
    io.add_object('arm1', [Contactor('Mass1', contact_type='Face',
                                     contact_index=0),
                           Contactor('Arm1',
                                     contact_type='Face',
                                     contact_index=0,
                                     relative_translation=[0, r1+l1, 0],
                                     relative_orientation=[(1, 0, 0), pi/2])],
                  translation=[0, 0, r2 + gap + r2 + l2 + r1 + hgap],
                  orientation=((1, 0, 0), pi/2),
                  inertia=np.eye(3),
                  mass=m1)

    # second branch + second mass
    io.add_object('arm2', [Contactor('Mass2', instance_name='Mass2Contact',
                                     contact_type='Face', contact_index=0),
                           Contactor('Arm2',
                                     contact_type='Face',
                                     contact_index=0,
                                     relative_translation=[0, r2+l2, 0],
                                     relative_orientation=[(1, 0, 0), pi/2])],
                  translation=[0, 0, r2 + gap],
                  orientation=((1, 0, 0), pi/2),
                  velocity=[0, 20, 0, 0, 0, 0],
                  inertia=np.eye(3),
                  mass=m2)

    io.add_joint('joint1', 'arm1', 'arm2',
                 points=[[0, 0, -r1]],
                 axes=[[1, 0, 0]],
                 joint_class='PivotJointR',
                 absolute=False)

    io.add_joint('joint2', 'arm1',
                 points=[[0, r2 + gap + r2 + l2 + r1 + hgap + l1, 0]],
                 axes=[[1, 0, 0]],
                 joint_class='PivotJointR',
                 absolute=False)

    io.add_object('ball', [Contactor('Brick', instance_name='BallContact',
                                     contact_type='Face', contact_index=0)],
                  translation=[0, -2*r2, r2+gap+.5],
                  inertia=np.eye(3),
                  mass=m2)

    # the ground object made with the ground shape. As the mass is
    # not given, it is a static object only involved in contact
    # detection.
    io.add_object('ground', [Contactor('Ground',
                                       contact_type='Face', contact_index=5)],
                  translation=[0, 0, -.25])


    io.add_interaction('sphere-ball',
                       'arm2', 'Mass2Contact',
                       'ball', 'BallContact',
                       distance_calculator='occ',
                       offset1=0.01)

    io.add_interaction('ball-ground',
                       'ball', 'BallContact',
                       'ground', 'Ground-0',
                       distance_calculator='occ',
                       offset1=0.01)

# Run the simulation from the inputs previously defined and add
# results to the hdf5 file. The visualisation of the output may be done
# with the vview command.
with MechanicsHdf5Runner(mode='r+') as io:

    io.run()
