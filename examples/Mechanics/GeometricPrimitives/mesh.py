#!/usr/bin/env python

#
# Example of one object under gravity with one contactor and a ground
# using the Siconos proposed mechanics API
#

from siconos.mechanics.collision.tools import Contactor
from siconos.io.mechanics_io import Hdf5
import siconos.numerics as Numerics

# Creation of the hdf5 file for input/output
with Hdf5() as io:

    # Load a mesh.  The example mesh is a low-poly version of the
    # Stanford Bunny by Saf, license: Creative Commons - Attribution.
    # Taken from http://www.thingiverse.com/thing:466857
    io.addMeshFromFile('Bunny', "bunny.stl", scale=0.01,
                       insideMargin=0.0, outsideMargin=0.0)

    # Definition of the ground shape
    io.addPrimitiveShape('Ground', 'Box', (10, 10, 1.0))

    # Definition of a non smooth law. As no group ids are specified it
    # is between contactors of group id 0.
    io.addNewtonImpactFrictionNSL('contact', mu=0.1, e=0.7)

    # The sphere object made with an unique Contactor : the sphere shape.
    # As a mass is given, it is a dynamic system involved in contact
    # detection and in the simulation.  With no group id specified the
    # Contactor belongs to group 0
    io.addObject('bunny', [Contactor('Bunny')],
                 translation=[0, 0, 4],
                 orientation=[1, 0, 0, 0],
                 velocity=[0, 0, 0, 1, 0, 0],
                 mass=1)

    # the ground object made with the ground shape. As the mass is
    # not given, it is a static object only involved in contact
    # detection.
    io.addObject('ground', [Contactor('Ground')],
                 translation=[0, 0, -0.5])

# Run the simulation from the inputs previously defined and add
# results to the hdf5 file. The visualisation of the output may be done
# with the vview command.
with Hdf5(mode='r+') as io:

    # By default earth gravity is applied and the units are those
    # of the International System of Units.
    io.run(with_timer=False,
           t0=0,
           T=20,
           h=0.005,
           multipoints_iterations=True,
           theta=0.50001,
           Newton_max_iter=20,
           set_external_forces=None,
           solver=Numerics.SICONOS_FRICTION_3D_NSGS,
           itermax=100000,
           tolerance=1e-8,
           numerics_verbose=False,
           output_frequency=None)
