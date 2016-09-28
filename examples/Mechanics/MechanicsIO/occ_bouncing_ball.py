#!/usr/bin/env python

#
# Example of a bouncing ball with OpenCascade contactors
#

from siconos.mechanics.contact_detection.tools import Contactor
from siconos.mechanics import occ
from siconos.io.mechanics_io import Hdf5
import siconos.numerics as Numerics

from OCC.BRepPrimAPI import BRepPrimAPI_MakeBox, BRepPrimAPI_MakeSphere

sphere = BRepPrimAPI_MakeSphere(1.).Shape()
ground = BRepPrimAPI_MakeBox(100., 100., .5).Shape()

# Creation of the hdf5 file for input/output
with Hdf5() as io:

    io.addOccShape('Sphere', sphere)
    io.addOccShape('Ground', ground)

    io.addObject('sphere',
                 [Contactor('Sphere', contact_type='Face', contact_index=0)],
                 mass=1, translation=[50, 50, 10], velocity=[0, 0, 0, 0, 0, 0])

    io.addObject('ground',
                 [Contactor('Ground', contact_type='Face', contact_index=5)],
                 mass=0, translation=[0, 0, 0])

    io.addPermanentInteraction('sphere-ground',
                               'sphere', 'Sphere-0',
                               'ground', 'Ground-0')

    io.addNewtonImpactFrictionNSL('contact', mu=0.3, e=0.9)

# Run the simulation from the inputs previously defined and add
# results to the hdf5 file. The visualisation of the output may be done
# with the vview command.
with Hdf5(mode='r+') as io:

    # By default earth gravity is applied and the units are those
    # of the International System of Units.
    # Because of fixed collision margins used in the collision detection,
    # sizes of small objects may need to be expressed in cm or mm.

    io.run(with_timer=False,
           time_stepping=occ.OccTimeStepping,
           space_filter=occ.OccSpaceFilter,
           body_class=None,
           shape_class=None,
           face_class=None,
           edge_class=None,
           gravity_scale=1,
           t0=0,
           T=10,
           h=0.0005,
           multipoints_iterations=True,
           theta=0.50001,
           Newton_max_iter=20,
           set_external_forces=None,
           solver=Numerics.SICONOS_FRICTION_3D_NSGS,
           itermax=100000,
           tolerance=1e-8,
           numerics_verbose=False,
           output_frequency=None)
