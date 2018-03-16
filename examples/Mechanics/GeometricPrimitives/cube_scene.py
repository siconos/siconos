#!/usr/bin/env python
from __future__ import print_function

#
# Example of one object under gravity with one contactor and a ground
#

from siconos.mechanics.collision.tools import Contactor
from siconos.io.mechanics_hdf5 import MechanicsHdf5

# Creation of the hdf5 file for input/output
with MechanicsHdf5(io_filename='cube_scene.hdf5') as io:

    # Definition of a cube as a convex shape
    io.add_convex_shape('Cube', [
        (-1.0, 1.0, -1.0),
        (-1.0, -1.0, -1.0),
        (-1.0, -1.0, 1.0),
        (-1.0, 1.0, 1.0),
        (1.0, 1.0, 1.0),
        (1.0, 1.0, -1.0),
        (1.0, -1.0, -1.0),
        (1.0, -1.0, 1.0)])

    # Alternative to the previous convex shape definition.
    # io.add_primitive_shape('Cube1', 'Box', (2, 2, 2))

    # Definition of the ground shape
    io.add_primitive_shape('Ground', 'Box', (100, 100, .5))

    # Definition of a non smooth law. As no group ids are specified it
    # is between contactors of group id 0.
    io.add_Newton_impact_friction_nsl('contact', mu=0.3)

    # The cube object made with an unique Contactor : the cube shape.
    # As a mass is given, it is a dynamic system involved in contact
    # detection and in the simulation.  With no group id specified the
    # Contactor belongs to group 0
    io.add_object('cube', [Contactor('Cube')], translation=[0, 0, 2],
                  velocity=[10, 0, 0, 1, 1, 1],
                  mass=1)

    # the ground object made with the ground shape. As the mass is
    # not given, it is a static object only involved in contact
    # detection.
    io.add_object('ground', [Contactor('Ground')],
                  translation=[0, 0, 0])

    print('Wrote', io._io_filename)
    print('Now run cube_simulation.py.')
