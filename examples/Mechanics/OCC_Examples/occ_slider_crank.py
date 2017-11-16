#!/usr/bin/env python

from siconos.mechanics.collision.tools import Volume, Contactor, Shape
from siconos.io.mechanics_io import Hdf5
import siconos.io.mechanics_io

siconos.io.mechanics_io.set_implementation('original')
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

    io.addPluginSource('plugin', 'SliderCrankPlugin/SliderCrankPlugin.cpp')

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

    io.addShapeDataFromFile('Artefact',
                            '../Mechanisms/SliderCrank/CAD/artefact2.step')

#    io.addObject('Artefact', [Shape(shape_data='Artefact',
#                                    instance_name='artefact')],
#                 translation=[0., 0., 0.])

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

    io.addObject('slider', [
        Shape(shape_data='Slider',
              instance_name='cslid',
              relative_translation=[-a, 0., 0.]),
        Contactor(
            instance_name='Contact_b_f1',
            shape_data='Contact_b_cyl',
            contact_type='Face',
            contact_index=1,
            relative_translation=[-a, 0., 0.]),
        Contactor(
            instance_name='Contact_h_f1',
            shape_data='Contact_h_cyl',
            contact_type='Face',
            contact_index=1,
            relative_translation=[-a, 0., 0.]),
        Contactor(
            instance_name='Contact_b_f0',
            shape_data='Contact_b_cyl',
            contact_type='Face',
            contact_index=0,
            relative_translation=[-a, 0., 0.]),
        Contactor(
            instance_name='Contact_h_f0',
            shape_data='Contact_h_cyl',
            contact_type='Face',
            contact_index=0,
            relative_translation=[-a, 0., 0.]),
        Contactor(
            instance_name='Contact_b_e1',
            shape_data='Contact_b_cyl',
            contact_type='Edge',
            contact_index=1,
            relative_translation=[-a, 0., 0.]),
        Contactor(
            instance_name='Contact_h_e1',
            shape_data='Contact_h_cyl',
            contact_type='Edge',
            contact_index=1,
            relative_translation=[-a, 0., 0.]),
        Contactor(
            instance_name='Contact_b_e0',
            shape_data='Contact_b_cyl',
            contact_type='Edge',
            contact_index=0,
            relative_translation=[-a, 0., 0.]),
        Contactor(
            instance_name='Contact_h_e0',
            shape_data='Contact_h_cyl',
            contact_type='Edge',
            contact_index=0,
            relative_translation=[-a, 0., 0.])],
                 translation=[l1 + l2 + a, 0., 0.],
                 velocity=[-0., 0., 0., 0., w30, 0.],
                 mass=0.076,
                 inertia=[2.7e-6, 1., 1.])

    # a static object (mass=0)
    io.addObject('chamber', [Contactor(
        instance_name='Chamber_contact0',
        shape_data='Chamber',
        contact_type='Face',
        contact_index=0,
        relative_translation=[0, 0, 0]),
                            Contactor(
        instance_name='Chamber_contact1',
        shape_data='Chamber',
        contact_type='Face',
        contact_index=1,
        relative_translation=[0, 0, 0])],
                translation=[0, 0, 0])

    io.addJoint('joint1',  'part1',
                points=[[0., 0., 0.]],
                axes=[[0., 1., 0.]],
                joint_class='PivotJointR')

    io.addJoint('joint2', 'part2', 'slider',
                pivot_point=[0.5*l2, 0., 0.],
                axis=[0., 1., 0],
                joint_class='PivotJointR')

    io.addJoint('joint3', 'part1', 'part2',
                pivot_point=[l1, 0., 0.],
                axis=[0., 1., 0.],
                joint_class='PivotJointR')

    io.addInteraction('contact10',
                      body1_name='slider', contactor1_name='Contact_b_f0',
                      body2_name='chamber', contactor2_name='Chamber_contact0',
                      distance_calculator='cadmbtb',
                      offset=0.024)

    io.addInteraction('contact11',
                      body1_name='slider', contactor1_name='Contact_b_f0',
                      body2_name='chamber', contactor2_name='Chamber_contact1',
                      distance_calculator='cadmbtb',
                      offset=0.024)

    io.addInteraction('contact20',
                      body1_name='slider', contactor1_name='Contact_h_f0',
                      body2_name='chamber', contactor2_name='Chamber_contact0',
                      distance_calculator='cadmbtb',
                      offset=0.024)

    io.addInteraction('contact21',
                      body1_name='slider', contactor1_name='Contact_h_f0',
                      body2_name='chamber', contactor2_name='Chamber_contact1',
                      distance_calculator='cadmbtb',
                      offset=0.024)

    io.addInteraction('contact30',
                      body1_name='slider', contactor1_name='Contact_b_f1',
                      body2_name='chamber', contactor2_name='Chamber_contact0',
                      distance_calculator='cadmbtb',
                      offset=0.024)

    io.addInteraction('contact31',
                      body1_name='slider', contactor1_name='Contact_b_f1',
                      body2_name='chamber', contactor2_name='Chamber_contact1',
                      distance_calculator='cadmbtb',
                      offset=0.024)

    io.addInteraction('contact40',
                      body1_name='slider', contactor1_name='Contact_h_f1',
                      body2_name='chamber', contactor2_name='Chamber_contact0',
                      distance_calculator='cadmbtb',
                      offset=0.024)

    io.addInteraction('contact41',
                      body1_name='slider', contactor1_name='Contact_h_f1',
                      body2_name='chamber', contactor2_name='Chamber_contact1',
                      distance_calculator='cadmbtb',
                      offset=0.024)

    io.addInteraction('contact50',
                      body1_name='slider', contactor1_name='Contact_b_e0',
                      body2_name='chamber', contactor2_name='Chamber_contact0',
                      distance_calculator='cadmbtb',
                      offset=0.024)

    io.addInteraction('contact51',
                      body1_name='slider', contactor1_name='Contact_b_e0',
                      body2_name='chamber', contactor2_name='Chamber_contact1',
                      distance_calculator='cadmbtb',
                      offset=0.024)

    io.addInteraction('contact60',
                      body1_name='slider', contactor1_name='Contact_h_e0',
                      body2_name='chamber', contactor2_name='Chamber_contact0',
                      distance_calculator='cadmbtb',
                      offset=0.024)

    io.addInteraction('contact61',
                      body1_name='slider', contactor1_name='Contact_h_e0',
                      body2_name='chamber', contactor2_name='Chamber_contact1',
                      distance_calculator='cadmbtb',
                      offset=0.024)

    io.addInteraction('contact70',
                      body1_name='slider', contactor1_name='Contact_b_e1',
                      body2_name='chamber', contactor2_name='Chamber_contact0',
                      distance_calculator='cadmbtb',
                      offset=0.024)

    io.addInteraction('contact71',
                      body1_name='slider', contactor1_name='Contact_b_e1',
                      body2_name='chamber', contactor2_name='Chamber_contact1',
                      distance_calculator='cadmbtb',
                      offset=0.024)

    io.addInteraction('contact80',
                      body1_name='slider', contactor1_name='Contact_h_e1',
                      body2_name='chamber', contactor2_name='Chamber_contact0',
                      distance_calculator='cadmbtb',
                      offset=0.024)

    io.addInteraction('contact81',
                      body1_name='slider', contactor1_name='Contact_h_e1',
                      body2_name='chamber', contactor2_name='Chamber_contact1',
                      distance_calculator='cadmbtb',
                      offset=0.024)

    io.addExternalFunction('f1', 'part1', 'setComputeFExtFunction',
                           'SliderCrankPlugin', 'externalForcesB1')

    io.addExternalFunction('f2', 'part2', 'setComputeFExtFunction',
                           'SliderCrankPlugin', 'externalForcesB2')

    io.addExternalFunction('f3', 'slider', 'setComputeFExtFunction',
                           'SliderCrankPlugin', 'externalForcesS')

    io.addExternalBCFunction('fbc', 'part1', [4],
                             'SliderCrankPlugin', 'prescribedvelocityB1')

    io.addNewtonImpactFrictionNSL('contact', mu=0.3, e=0.4)

with Hdf5(mode='r+') as io:

    io.run(with_timer=True,
           t0=0,
           T=10,
           h=0.0005,
           Newton_max_iter=5)
