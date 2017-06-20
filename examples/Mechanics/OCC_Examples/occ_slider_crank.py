#!/usr/bin/env python

from siconos.mechanics.collision.tools import Volume, Contactor
from siconos.io.mechanics_io import Hdf5
import siconos.io.mechanics_io

from OCC.BRepPrimAPI import BRepPrimAPI_MakeBox, BRepPrimAPI_MakeSphere
from OCC.gp import gp_Pnt

siconos.io.mechanics_io.set_implementation('original')
siconos.io.mechanics_io.set_backend('occ')

sphere = BRepPrimAPI_MakeSphere(0.025).Shape()

l1 = 0.153    # crank length
l2 = 0.306    # connecting rod length
a = 0.05      # half length of the slider
b = 0.025     # half height of the slider
c = 0.001     # clearance between slider and guide
w10 = -150.   # initial angular speed for the crank
w20 = 75.     # initial angular speed for the connecting rod
w30 = 0.      # initial angular speed for the slider

with Hdf5() as io:

    io.addOccShape('Sphere', sphere)
    
    io.addShapeDataFromFile('body1',
                            '../Mechanisms/SliderCrank/CAD/body1.step')
    io.addShapeDataFromFile('body2',
                            '../Mechanisms/SliderCrank/CAD/body2.step')
    io.addShapeDataFromFile('Slider',
                            '../Mechanisms/SliderCrank/CAD/Slider.step')

    io.addShapeDataFromFile('Contact_b_cyl',
                            '../Mechanisms/SliderCrank/CAD/contact_b_cyl.step')

    io.addShapeDataFromFile('Contact_h_cyl',
                            '../Mechanisms/SliderCrank/CAD/contact_h_cyl.step')

    io.addShapeDataFromFile('RingBody',
                            '../Mechanisms/SliderCrank/CAD/RingBody1.stp')

    io.addShapeDataFromFile('Chamber',
                            '../Mechanisms/SliderCrank/CAD/chamber.step')

    io.addShapeDataFromFile('AxisBody',
                            '../Mechanisms/SliderCrank/CAD/AxisBody2.stp')

    io.addObject('part1', [Volume(shape_data='body1',
                                  instance_name='Body1',
                                  relative_translation=[-0.5*l1, 0., 0.],
                                  relative_orientation=[(0, 1, 0), 0. ])],
                 translation=[0.5*l1, 0., 0.],
                 velocity=[0., 0., -0.5 * w10 * l1, 0., w10, 0.],
                 mass=0.038,
                 inertia=[7.4e-5, 1, 1.])

    io.addObject('part2', [Volume(shape_data='body2',
                                  instance_name='Body2',
                                  relative_translation=[-0.5 * l2, 0., 0.])],
                 translation=[l1 + 0.5*l2, 0., 0.],
                 orientation=[0., 0., 1., 0.],
                 velocity=[0., 0., -0.5 * w10 * l1, 0., w20, 0.],
                 mass=0.038,
                 inertia=[5.9e-4, 1., 1.])

    io.addObject('slider', [Contactor(shape_data='Sphere',
                                      instance_name='cslid',
                                      contact_type='Face',
                                      contact_index=0,
                                      relative_translation=[-a, 0., 0.])],
                            # Contactor(
                            #     instance_name='Contact_b',
                            #     shape_data='Contact_b_cyl',
                            #     contact_type='Edge',
                            #     contact_index=0,
                            #     relative_translation=[-a, 0., 0.]),
                            # Contactor(
                            #     instance_name='Contact_h',
                            #     shape_data='Contact_h_cyl',
                            #     contact_type='Edge',
                            #     contact_index=0,
                            #     relative_translation=[-a, 0., 0.])],
                 translation=[l1 + l2 + a, 0., 0.],
                 velocity=[0., 0., 0., 0., w30, 0.],
                 mass=0.076,
                 inertia=[2.7e-6, 1., 1.])

    # a static object (mass=0)
    io.addObject('chamber', [Contactor(
        instance_name='Chamber_contact',
        shape_data='Chamber',
        contact_type='Face',
        contact_index=0,
        relative_translation=[0, 0, 0])],
                 translation=[0, 0, 0])

    io.addJoint('joint1',  'part1',
                pivot_point=[0., 0., 0.],
                axis=[0., 1., 0.],
                joint_class='PivotJointR')

    io.addJoint('joint2', 'part2', 'slider',
                pivot_point=[-0.5*l2, 0., 0.],
                axis=[0., 1., 0],
                joint_class='PivotJointR')

    io.addJoint('joint3', 'part1', 'part2',
                pivot_point=[0.5*l1, 0., 0.],
                axis=[0., 1., 0.],
                joint_class='PivotJointR')

    io.addInteraction('contact1',
                      body1_name='slider', contactor1_name='cslid',
                      body2_name='chamber', contactor2_name='Chamber_contact',
                      distance_calculator='cadmbtb',
                      offset=0.0024)

    io.addNewtonImpactFrictionNSL('contact', mu=0.3, e=0.9)
    
with Hdf5(mode='r+') as io:

    io.run(with_timer=True,
           t0=0,
           T=1,
           h=0.0005,
           Newton_max_iter=1,
           set_external_forces=lambda x: None)
