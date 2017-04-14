#!/usr/bin/env python

#
# A simple walking mechanism (in progress, does not work!)
#

from siconos.mechanics.collision.tools import Contactor
from siconos.io.mechanics_io import Hdf5
import siconos.numerics as Numerics
import math, numpy

pi = math.pi

motor_id = None

# Creation of the hdf5 file for input/output
with Hdf5() as io:

    # Only touch the ground, ignore contacts between links of the robot
    io.addNewtonImpactFrictionNSL('contact',
                                  collision_group1 = 0,
                                  collision_group2 = 1,
                                  mu=0.3, e=0.0)

    # Definition of the ground shape
    io.addPrimitiveShape('Ground', 'Box', (100, 100, 2))
    io.addObject('ground', [Contactor('Ground', collision_group = 1)],
                 translation=[0, 0, -1])

    # Define shape of a bar-shaped link
    io.addPrimitiveShape('Bar1', 'Box', (10, 1, 1))
    io.addPrimitiveShape('Bar2', 'Box', (12.67, 1, 1))

    ## Core
    bar1 = \
    io.addObject('bar1', [Contactor('Bar1'),
                          Contactor('Bar1',
                                    relative_translation=[5.5, 0, 0],
                                    relative_orientation=[(0,0,1), pi/2]),
                          Contactor('Bar1',
                                    relative_translation=[-5.5, 0, 0],
                                    relative_orientation=[(0,0,1), pi/2])],
                 translation=[0, 0, 11],
                 mass=1)

    ## Legs
    io.addObject('bar2', [Contactor('Bar1')], translation=[5.5, 5.5, 6],
                 orientation=[(0,1,0), -pi/2], mass=1)

    io.addObject('bar3', [Contactor('Bar1')], translation=[5.5, -5.5, 6],
                 orientation=[(0,1,0), -pi/2], mass=1)

    io.addObject('bar4', [Contactor('Bar1')], translation=[-5.5, 5.5, 6],
                 orientation=[(0,1,0), -pi/2], mass=1)

    io.addObject('bar5', [Contactor('Bar1')], translation=[-5.5, -5.5, 6],
                 orientation=[(0,1,0), -pi/2], mass=1)

    io.addJoint('joint1', 'bar2', 'bar1', [4.5, 0, 0], [0, 1, 0], 'PivotJointR')
    io.addJoint('joint2', 'bar3', 'bar1', [4.5, 0, 0], [0, 1, 0], 'PivotJointR')
    io.addJoint('joint3', 'bar4', 'bar1', [4.5, 0, 0], [0, 1, 0], 'PivotJointR')
    io.addJoint('joint4', 'bar5', 'bar1', [4.5, 0, 0], [0, 1, 0], 'PivotJointR')

    ## Stabilizing leg links to body

    # See mech.spbstu.ru/Dzenushko_Dainis:_Walking_mechanisms_survey
    # for a survey of better walking mechanisms

    io.addObject('bar6', [Contactor('Bar2')], translation=[0,  6.5, 7],
                 orientation=[(0,1,0), -pi/8], mass=1)
    io.addObject('bar7', [Contactor('Bar2')], translation=[0, -6.5, 7],
                 orientation=[(0,1,0), -pi/8], mass=1)
    io.addObject('bar8', [Contactor('Bar2')], translation=[0,  6.5, 7],
                 orientation=[(0,1,0), -pi/8], mass=1)
    io.addObject('bar9', [Contactor('Bar2')], translation=[0, -6.5, 7],
                 orientation=[(0,1,0), -pi/8], mass=1)
    io.addJoint('joint5', 'bar6', 'bar2', [ 4.5,  0, 0], [0, 1, 0], 'PivotJointR')
    io.addJoint('joint6', 'bar7', 'bar3', [ 4.5,  0, 0], [0, 1, 0], 'PivotJointR')
    io.addJoint('joint7', 'bar8', 'bar4', [-4.5,  0, 0], [0, 1, 0], 'PivotJointR')
    io.addJoint('joint8', 'bar9', 'bar5', [-4.5,  0, 0], [0, 1, 0], 'PivotJointR')

    io.addJoint('joint9', 'bar6', 'bar8', [1, 0, 0], [], 'PrismaticJointR')
    io.addJoint('joint10','bar7', 'bar9', [1, 0, 0], [], 'PrismaticJointR')

    io.addJoint('joint11','bar2', 'bar3', [0, 0, 0], [0, 1, 0], 'PivotJointR')
    io.addJoint('joint12','bar4', 'bar5', [0, 0, 0], [0, 1, 0], 'PivotJointR')

    # Harmonic oscillator on Y-axis angular velocity = a+b*cos(omega*time+phi))
    freq = 2
    amp = 0.2
    io.addBoundaryCondition('vibration', 'bar1',
                            indices=[4], # Y-angular axis is index 4 into ds->v
                            bc_class='HarmonicBC',
                            a =     [     0.0 ],
                            b =     [     amp ],
                            omega = [ pi*freq ],
                            phi =   [     0.0 ])

def my_forces(body):
    g = 9.81
    weight = numpy.array([0, 0, - body.scalarMass() * g, 0, 0, 0])
    twist = numpy.array([0,0,0,0,300,0])
    push = numpy.array([0,0,0,0,0,0])
    force = weight
    if body.number() == bar1:
        force = weight + push # + twist
    body.setFExtPtr(force)

# Run the simulation from the inputs previously defined and add
# results to the hdf5 file. The visualisation of the output may be done
# with the vview command.
with Hdf5(mode='r+') as io:

    io.run(with_timer=False,
            time_stepping=None,
            space_filter=None,
            body_class=None,
            shape_class=None,
            face_class=None,
            edge_class=None,
            t0=0,
            T=30,
            h=0.005,
            multipoints_iterations=True,
            theta=0.50001,
            Newton_max_iter=1,
            set_external_forces=my_forces,
            solver=Numerics.SICONOS_FRICTION_3D_NSGS,
            itermax=1000,
            tolerance=1e-4,
            numerics_verbose=False,
            output_frequency=None)
