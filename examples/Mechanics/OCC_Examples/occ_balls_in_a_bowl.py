#!/usr/bin/env python

#
# Tree balls in a spinning bowl, with gravity and a ground
# OpenCascade contactors
#
# specification of center of mass & moments of inertia

#

from siconos.mechanics.collision.tools import Contactor
from siconos.io.mechanics_run import MechanicsHdf5Runner
from siconos import numerics

# for osi specification:
# from siconos import kernel

from OCC.BRepAlgoAPI import BRepAlgoAPI_Cut
from OCC.BRepPrimAPI import BRepPrimAPI_MakeBox, BRepPrimAPI_MakeSphere
from OCC.gp import gp_Pnt, gp_Ax1, gp_Dir

from OCC.GProp import GProp_GProps
from OCC.BRepGProp import brepgprop_VolumeProperties

from math import pi

# original implementation with occ backend
import siconos.io.mechanics_run
siconos.io.mechanics_run.set_backend('occ')

# ball shape
ball = BRepPrimAPI_MakeSphere(.15).Shape()

ball_props = GProp_GProps()
brepgprop_VolumeProperties(ball, ball_props)

ball_mass = ball_props.Mass()
ball_com = ball_props.CentreOfMass()
ball_inertia = ball_props.MatrixOfInertia()

ball_I1 = ball_props.MomentOfInertia(gp_Ax1(ball_com, gp_Dir(1, 0, 0)))
ball_I2 = ball_props.MomentOfInertia(gp_Ax1(ball_com, gp_Dir(0, 1, 0)))
ball_I3 = ball_props.MomentOfInertia(gp_Ax1(ball_com, gp_Dir(0, 0, 1)))

print('ball mass:', ball_mass)
print('ball center of mass:', (ball_com.Coord(1),
                               ball_com.Coord(2),
                               ball_com.Coord(3)))
print('ball moment of inertia:', (ball_I1, ball_I2, ball_I3))

# the ground
ground = BRepPrimAPI_MakeBox(gp_Pnt(-20, -20, 0), 40., 40., .5).Shape()

# bowl shape
hspherei = BRepPrimAPI_MakeSphere(.9, pi).Shape()
hsphereo = BRepPrimAPI_MakeSphere(1., pi).Shape()
bowl = BRepAlgoAPI_Cut(hsphereo, hspherei).Shape()

bowl_props = GProp_GProps()
brepgprop_VolumeProperties(bowl, bowl_props)

bowl_mass = bowl_props.Mass()
bowl_com = bowl_props.CentreOfMass()
bowl_inertia = bowl_props.MatrixOfInertia()

bowl_I1 = bowl_props.MomentOfInertia(gp_Ax1(bowl_com, gp_Dir(1, 0, 0)))
bowl_I2 = bowl_props.MomentOfInertia(gp_Ax1(bowl_com, gp_Dir(0, 1, 0)))
bowl_I3 = bowl_props.MomentOfInertia(gp_Ax1(bowl_com, gp_Dir(0, 0, 1)))

print('bowl mass:', bowl_mass)
print('bowl center of mass:', (bowl_com.Coord(1),
                               bowl_com.Coord(2),
                               bowl_com.Coord(3)))
print('bowl moment of inertia:', (bowl_I1, bowl_I2, bowl_I3))

# Creation of the hdf5 file for input/output
with MechanicsHdf5Runner() as io:

    io.add_occ_shape('Contact', bowl)

    io.add_occ_shape('Ground', ground)

    io.add_occ_shape('Ball', ball)

    io.add_object('bowl',
                  [Contactor('Contact',
                             contact_type='Face',
                             contact_index=0,
                  ),
                   Contactor('Contact',
                             contact_type='Face',
                             contact_index=3,
                   ),
                   Contactor('Contact',
                             contact_type='Edge',
                             contact_index=0,
                   )],
                  mass=bowl_mass,

                  orientation=([1, 0, 0], -pi / 2),
                  translation=[0, 0, 2],
                  velocity=[0, 0, 0, 0, 2, 0],
                  center_of_mass=[bowl_com.Coord(1),
                                  bowl_com.Coord(2),
                                  bowl_com.Coord(3)],
                  inertia=[bowl_I1, bowl_I2, bowl_I3])

    #
    # balls
    #

    io.add_object('ball1',
                  [Contactor('Ball',
                             instance_name='Ball1',
                             contact_type='Face',
                             contact_index=0)],
                  translation=[0, .3, 2],
                  mass=.1,
                  inertia=[ball_I1, ball_I2, ball_I3])

    io.add_object('ball2',
                  [Contactor('Ball',
                             instance_name='Ball2',
                             contact_type='Face',
                             contact_index=0)],
                  translation=[0, 0, 2], mass=.1,
                  inertia=[ball_I1, ball_I2, ball_I3])

    io.add_object('ball3',
                  [Contactor('Ball',
                             instance_name='Ball3',
                             contact_type='Face',
                             contact_index=0)],
                  translation=[0, -.3, 2],
                  mass=.1,
                  inertia=[ball_I1, ball_I2, ball_I3])

    #
    # ground, static object (mass=0)
    #

    io.add_object('ground',
                  [Contactor('Ground',
                             contact_type='Face',
                             contact_index=5)],
                  mass=0,
                  translation=[0, 0, 0])

    #
    # interactions, order ball -> bowl is important
    # ball -> ground if some balls are ejected

    io.add_interaction('bowl-ground',
                       'bowl', 'Contact-0',
                       'ground', 'Ground-0',
                       distance_calculator='cadmbtb',
                       offset1=0.01)

    io.add_interaction('bowl-ball1',
                       'ball1', 'Ball1',
                       'bowl', 'Contact-1',
                       distance_calculator='cadmbtb',
                       offset1=0.05)

    io.add_interaction('bowl-ball2',
                       'ball2', 'Ball2',
                       'bowl', 'Contact-1',
                       distance_calculator='cadmbtb',
                       offset1=0.05)

    io.add_interaction('bowl-ball3',
                       'ball3', 'Ball3',
                       'bowl', 'Contact-1',
                       distance_calculator='cadmbtb',
                       offset1=0.05)

    io.add_interaction('ball1-ball2',
                       'ball1', 'Ball1',
                       'ball2', 'Ball2',
                       distance_calculator='cadmbtb',
                       offset1=0.05)

    io.add_interaction('ball1-ball3',
                       'ball1', 'Ball1',
                       'ball3', 'Ball3',
                       distance_calculator='cadmbtb',
                       offset1=0.05)

    io.add_interaction('ball2-ball3',
                       'ball2', 'Ball2',
                       'ball3', 'Ball3',
                       distance_calculator='cadmbtb',
                       offset1=0.05)

    io.add_Newton_impact_friction_nsl('contact', mu=0.3, e=0.)

# Run the simulation from the inputs previously defined and add
# results to the hdf5 file. The visualisation of the output may be done
# with the vview command.
with MechanicsHdf5Runner(mode='r+') as io:

    io.run(with_timer=False,
           gravity_scale=1,
           t0=0,
           T=10,
           h=0.0005,
           theta=0.50001,
           Newton_max_iter=20,
           set_external_forces=None,
           solver=numerics.SICONOS_FRICTION_3D_NSGS,
           itermax=100000,
           tolerance=1e-7,
           numerics_verbose=False,
           output_frequency=None
           # osi=kernel.MoreauJeanCombinedProjectionOSI
           )
