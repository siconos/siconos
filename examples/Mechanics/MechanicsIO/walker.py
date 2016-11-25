#!/usr/bin/env python

#
# A simple walking mechanism (in progress, does not work!)
#

from siconos.mechanics.collision.tools import Contactor
from siconos.io.mechanics_io import Hdf5
import siconos.numerics as Numerics
import math

motor_id = None

# Creation of the hdf5 file for input/output
with Hdf5() as io:
    io.addNewtonImpactFrictionNSL('contact', mu=0.3)

    # Definition of the ground shape
    io.addPrimitiveShape('Ground', 'Box', (100, 100, 2))
    io.addObject('ground', [Contactor('Ground')], translation=[0, 0, -1])

    # Define shape of a bar-shaped link
    io.addPrimitiveShape('Bar1', 'Box', (10, 1, 1))

    motor_id = io.addObject('bar1', [Contactor('Bar1')],
                            translation=[0, 0, 20], mass=1)

    io.addObject('bar2', [Contactor('Bar1')], translation=[5.5, 0, 20],
                 orientation=[(0,0,1),math.pi/2], mass=1)

    io.addObject('bar3', [Contactor('Bar1')], translation=[5.5, 5.5, 15],
                 orientation=[(0,1,0),math.pi/2], mass=1)

    io.addObject('bar4', [Contactor('Bar1')], translation=[5.5, -5.5, 15],
                 orientation=[(0,1,0),math.pi/2], mass=1)

    io.addJoint('joint2', 'bar1', 'bar2', [5.5,    0, 20], [0, 0, 1], 'PivotJointR')
    io.addJoint('joint3', 'bar2', 'bar3', [5.5,  5.5, 20], [1, 0, 0], 'PivotJointR')
    io.addJoint('joint4', 'bar2', 'bar4', [5.5, -5.5, 20], [1, 0, 0], 'PivotJointR')

def apply_force(body):
    if body.number() == motor_id:
        body.setFExtPtr([0,0,0,10,0,0])

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
            T=10,
            h=0.01,
            multipoints_iterations=True,
            theta=0.50001,
            Newton_max_iter=1,
            set_external_forces=apply_force,
            solver=Numerics.SICONOS_FRICTION_3D_NSGS,
            itermax=1000,
            tolerance=1e-4,
            numerics_verbose=False,
            output_frequency=None)
