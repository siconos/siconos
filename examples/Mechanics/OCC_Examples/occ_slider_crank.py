#!/usr/bin/env python

from siconos.mechanics.collision.tools import Volume, Contactor, Shape
from siconos.io.mechanics_io import Hdf5
import siconos.io.mechanics_io

siconos.io.mechanics_io.set_backend('occ')

l1 = 0.153    # crank length
l2 = 0.306    # connecting rod length
a = 0.05      # half length of the slider
b = 0.025     # half height of the slider
c = 0.001     # clearance between slider and guide
w10 = -150.   # initial angular speed for the crank
w20 = 75.     # initial angular speed for the connecting rod
w30 = 0.      # initial angular speed for the slider

with Hdf5() as io:

    io.add_plugin_source('plugin', 'SliderCrankPlugin/SliderCrankPlugin.cpp')

    io.add_shape_data_from_file('body1',
                                '../Mechanisms/SliderCrank/CAD/body1.step')
    io.add_shape_data_from_file('body2',
                                '../Mechanisms/SliderCrank/CAD/body2.step')
    io.add_shape_data_from_file('Slider',
                                '../Mechanisms/SliderCrank/CAD/Slider.step')

    io.add_shape_data_from_file(
        'Contact_b_cyl',
        '../Mechanisms/SliderCrank/CAD/contact_b_cyl.step')

    io.add_shape_data_from_file(
        'Contact_h_cyl',
        '../Mechanisms/SliderCrank/CAD/contact_h_cyl.step')
    
    io.add_shape_data_from_file('RingBody',
                                '../Mechanisms/SliderCrank/CAD/RingBody1.stp')
    
    io.add_shape_data_from_file('Chamber',
                                '../Mechanisms/SliderCrank/CAD/chamber.step')

    io.add_shape_data_from_file('AxisBody',
                                '../Mechanisms/SliderCrank/CAD/AxisBody2.stp')

    io.add_shape_data_from_file('Artefact',
                                '../Mechanisms/SliderCrank/CAD/artefact2.step')

    io.add_object('Artefact', [Shape(shape_name='Artefact',
                                    instance_name='artefact')],
                 translation=[0., 0., 0.])

    io.add_object('part1', [Volume(shape_name='body1',
                                  instance_name='Body1',
                                  relative_translation=[-0.5*l1, 0., 0.],
                                  relative_orientation=[(0, 1, 0), 0.])],
                 translation=[0.5*l1, 0., 0.],
                 velocity=[0., 0., -0.5 * w10 * l1, 0., w10, 0.],
                 mass=0.038,
                 inertia=[7.4e-5, 1, 1.])

    io.add_object('part2', [Volume(shape_name='body2',
                                  instance_name='Body2',
                                  relative_translation=[-0.5 * l2, 0., 0.])],
                 translation=[l1 + 0.5*l2, 0., 0.],
                 orientation=[0., 0., 1., 0.],
                 velocity=[0., 0., -0.5 * w10 * l1, 0., w20, 0.],
                 mass=0.038,
                 inertia=[5.9e-4, 1., 1.])

    io.add_object('slider', [
        Shape(shape_name='Slider',
              instance_name='cslid',
              relative_translation=[-a, 0., 0.]),
        Contactor(
            instance_name='Contact_b_f1',
            shape_name='Contact_b_cyl',
            contact_type='Face',
            contact_index=1,
            relative_translation=[-a, 0., 0.]),
        Contactor(
            instance_name='Contact_h_f1',
            shape_name='Contact_h_cyl',
            contact_type='Face',
            contact_index=1,
            relative_translation=[-a, 0., 0.]),
        Contactor(
            instance_name='Contact_b_f0',
            shape_name='Contact_b_cyl',
            contact_type='Face',
            contact_index=0,
            relative_translation=[-a, 0., 0.]),
        Contactor(
            instance_name='Contact_h_f0',
            shape_name='Contact_h_cyl',
            contact_type='Face',
            contact_index=0,
            relative_translation=[-a, 0., 0.]),
        Contactor(
            instance_name='Contact_b_e1',
            shape_name='Contact_b_cyl',
            contact_type='Edge',
            contact_index=1,
            relative_translation=[-a, 0., 0.]),
        Contactor(
            instance_name='Contact_h_e1',
            shape_name='Contact_h_cyl',
            contact_type='Edge',
            contact_index=1,
            relative_translation=[-a, 0., 0.]),
        Contactor(
            instance_name='Contact_b_e0',
            shape_name='Contact_b_cyl',
            contact_type='Edge',
            contact_index=0,
            relative_translation=[-a, 0., 0.]),
        Contactor(
            instance_name='Contact_h_e0',
            shape_name='Contact_h_cyl',
            contact_type='Edge',
            contact_index=0,
            relative_translation=[-a, 0., 0.])],
                 translation=[l1 + l2 + a, 0., 0.],
                 velocity=[-0., 0., 0., 0., w30, 0.],
                 mass=0.076,
                 inertia=[2.7e-6, 1., 1.])

    # a static object (mass=0)
    io.add_object('chamber', [Contactor(
        instance_name='Chamber_contact0',
        shape_name='Chamber',
        contact_type='Face',
        contact_index=0,
        relative_translation=[0, 0, 0]),
                            Contactor(
        instance_name='Chamber_contact1',
        shape_name='Chamber',
        contact_type='Face',
        contact_index=1,
        relative_translation=[0, 0, 0])],
                translation=[0, 0, 0])

    io.add_joint('joint1',  'part1',
                 points=[[0., 0., 0.]],
                 axes=[[0., 1., 0.]],
                 joint_class='PivotJointR',
                 absolute=True)

    io.add_joint('joint2', 'part2', 'slider',
                 points=[[l1+l2, 0., 0.]],
                 axes=[[0., 1., 0]],
                 joint_class='PivotJointR',
                 absolute=True)

    io.add_joint('joint3', 'part1', 'part2',
                 points=[[l1, 0., 0.]],
                 axes=[[0., 1., 0.]],
                 joint_class='PivotJointR',
                 absolute=True)

    io.add_interaction('contact10',
                      body1_name='slider', contactor1_name='Contact_b_f0',
                      body2_name='chamber', contactor2_name='Chamber_contact0',
                      distance_calculator='cadmbtb',
                      offset=0.024)

    io.add_interaction('contact11',
                      body1_name='slider', contactor1_name='Contact_b_f0',
                      body2_name='chamber', contactor2_name='Chamber_contact1',
                      distance_calculator='cadmbtb',
                      offset=0.024)

    io.add_interaction('contact20',
                      body1_name='slider', contactor1_name='Contact_h_f0',
                      body2_name='chamber', contactor2_name='Chamber_contact0',
                      distance_calculator='cadmbtb',
                      offset=0.024)

    io.add_interaction('contact21',
                      body1_name='slider', contactor1_name='Contact_h_f0',
                      body2_name='chamber', contactor2_name='Chamber_contact1',
                      distance_calculator='cadmbtb',
                      offset=0.024)

    io.add_interaction('contact30',
                      body1_name='slider', contactor1_name='Contact_b_f1',
                      body2_name='chamber', contactor2_name='Chamber_contact0',
                      distance_calculator='cadmbtb',
                      offset=0.024)

    io.add_interaction('contact31',
                      body1_name='slider', contactor1_name='Contact_b_f1',
                      body2_name='chamber', contactor2_name='Chamber_contact1',
                      distance_calculator='cadmbtb',
                      offset=0.024)

    io.add_interaction('contact40',
                      body1_name='slider', contactor1_name='Contact_h_f1',
                      body2_name='chamber', contactor2_name='Chamber_contact0',
                      distance_calculator='cadmbtb',
                      offset=0.024)

    io.add_interaction('contact41',
                      body1_name='slider', contactor1_name='Contact_h_f1',
                      body2_name='chamber', contactor2_name='Chamber_contact1',
                      distance_calculator='cadmbtb',
                      offset=0.024)

    io.add_interaction('contact50',
                      body1_name='slider', contactor1_name='Contact_b_e0',
                      body2_name='chamber', contactor2_name='Chamber_contact0',
                      distance_calculator='cadmbtb',
                      offset=0.024)

    io.add_interaction('contact51',
                      body1_name='slider', contactor1_name='Contact_b_e0',
                      body2_name='chamber', contactor2_name='Chamber_contact1',
                      distance_calculator='cadmbtb',
                      offset=0.024)

    io.add_interaction('contact60',
                      body1_name='slider', contactor1_name='Contact_h_e0',
                      body2_name='chamber', contactor2_name='Chamber_contact0',
                      distance_calculator='cadmbtb',
                      offset=0.024)

    io.add_interaction('contact61',
                      body1_name='slider', contactor1_name='Contact_h_e0',
                      body2_name='chamber', contactor2_name='Chamber_contact1',
                      distance_calculator='cadmbtb',
                      offset=0.024)

    io.add_interaction('contact70',
                      body1_name='slider', contactor1_name='Contact_b_e1',
                      body2_name='chamber', contactor2_name='Chamber_contact0',
                      distance_calculator='cadmbtb',
                      offset=0.024)

    io.add_interaction('contact71',
                      body1_name='slider', contactor1_name='Contact_b_e1',
                      body2_name='chamber', contactor2_name='Chamber_contact1',
                      distance_calculator='cadmbtb',
                      offset=0.024)

    io.add_interaction('contact80',
                      body1_name='slider', contactor1_name='Contact_h_e1',
                      body2_name='chamber', contactor2_name='Chamber_contact0',
                      distance_calculator='cadmbtb',
                      offset=0.024)

    io.add_interaction('contact81',
                      body1_name='slider', contactor1_name='Contact_h_e1',
                      body2_name='chamber', contactor2_name='Chamber_contact1',
                      distance_calculator='cadmbtb',
                      offset=0.024)

    io.add_external_function('f1', 'part1', 'setComputeFExtFunction',
                           'SliderCrankPlugin', 'externalForcesB1')

    io.add_external_function('f2', 'part2', 'setComputeFExtFunction',
                           'SliderCrankPlugin', 'externalForcesB2')

    io.add_external_function('f3', 'slider', 'setComputeFExtFunction',
                           'SliderCrankPlugin', 'externalForcesS')

    io.add_external_bc_function('fbc', 'part1', [4],
                                'SliderCrankPlugin', 'prescribedvelocityB1')

    io.add_Newton_impact_friction_nsl('contact', mu=0.3, e=0.4)

with Hdf5(mode='r+') as io:

    io.run(with_timer=True,
           t0=0,
           T=1,
           h=0.0005,
           Newton_max_iter=5)
