#!/usr/bin/env python

# Various object types sliding, rolling, and sitting still.

from siconos.mechanics.collision.tools import Contactor
from siconos.io.mechanics_io import Hdf5
import siconos.numerics as Numerics
import siconos.io.mechanics_io
from siconos.mechanics.collision.convexhull import ConvexHull

options = siconos.mechanics.collision.bullet.SiconosBulletOptions()
options.worldScale = 1.0

# Control the number of perturbations applied to generate multipoint
# surface-surface contact manifolds.  Default is 5 and 5.
options.perturbationIterations = 5
options.minimumPointsPerturbationThreshold = 5

# Creation of the hdf5 file for input/output
with Hdf5() as io:

    # Definition of a tetrahedron as a convex shape.
    # Bottom purposely not even.
    import numpy
    pts = numpy.array([(-1.0, 1.0, -1.0),
                       (1.0, -1.0, -0.5),
                       (-1.0, -1.0, -0.7),
                       (0.0, 0.0, 1.0)])
    io.addConvexShape('Tetra', pts - pts.mean(0))

    io.addPrimitiveShape('Cyl', 'Cylinder', [1, 1])

    io.addPrimitiveShape('Cube', 'Box', [1, 1, 1])

    io.addPrimitiveShape('Ball', 'Sphere', [1])

    # Definition of the ground shape
    io.addPrimitiveShape('Ground', 'Box', (100, 100, 1))

    # Definition of a non smooth law. As no group ids are specified it
    # is between contactors of group id 0.
    io.addNewtonImpactFrictionNSL('contact', mu=0.3, e=0.0,
                                  collision_group1=1,
                                  collision_group2=2)

    # computation of inertia and volume
    ch = ConvexHull(pts)
    inertia,volume=ch.inertia(ch.centroid())

    # Copies of each object type, still and thrown horizontally.
    x = -20
    y = -40
    vel = 20
    spacing = 4

    io.addObject('tetra1', [Contactor('Tetra', collision_group=1)],
                 translation=[x, y, 2],
                 velocity=[0, 0, 0, 0, 0, 0],
                 mass=1, inertia=inertia)

    x += spacing
    io.addObject('tetra2', [Contactor('Tetra', collision_group=1)],
                 translation=[x, y, 2],
                 velocity=[0, vel, 0, 0, 0, 0],
                 mass=1, inertia=inertia)

    x += spacing
    io.addObject('cyl1', [Contactor('Cyl', collision_group=1)],
                 translation=[x, y, 2],
                 velocity=[0, 0, 0, 0, 0, 0],
                 mass=1)

    x += spacing
    io.addObject('cyl2', [Contactor('Cyl', collision_group=1)],
                 translation=[x, y, 2], orientation=[(1,1,1), numpy.pi/4],
                 velocity=[0, 0, 0, 0, 0, 1],
                 mass=1)

    x += spacing
    io.addObject('cyl3', [Contactor('Cyl', collision_group=1)],
                 translation=[x, y, 2],
                 velocity=[0, vel, 0, 0, 0, 0],
                 mass=1)

    x += spacing
    io.addObject('cyl4', [Contactor('Cyl', collision_group=1)],
                 translation=[x, y, 2], orientation=[(0,0,1), numpy.pi/2],
                 velocity=[0, vel, 0, 0, 0, 0],
                 mass=1)

    x += spacing
    io.addObject('ball1', [Contactor('Ball', collision_group=1)],
                 translation=[x, y, 2],
                 velocity=[0, 0, 0, 0, 0, 0],
                 mass=1)

    x += spacing
    io.addObject('ball2', [Contactor('Ball', collision_group=1)],
                 translation=[x, y, 2],
                 velocity=[0, vel, 0, 0, 0, 0],
                 mass=1)

    x += spacing
    io.addObject('box1', [Contactor('Cube', collision_group=1)],
                 translation=[x, y, 2],
                 velocity=[0, 0, 0, 0, 0, 0],
                 mass=1)

    x += spacing
    io.addObject('box2', [Contactor('Cube', collision_group=1)],
                 translation=[x, y, 2],
                 velocity=[0, vel, 0, 0, 0, 0],
                 mass=1)

    x += spacing
    stack_height = 3
    for i in range(stack_height):
        io.addObject('box%d'%(i+3), [Contactor('Cube', collision_group=1+(i%2))],
                     translation=[x, y, i+1],
                     velocity=[0, 0, 0, 0, 0, 0],
                     mass=1)

    # the ground object made with the ground shape. As the mass is
    # not given, it is a static object only involved in contact
    # detection.
    io.addObject('ground', [Contactor('Ground', collision_group=2)],
                 translation=[0, 0, 0])

# Run the simulation from the inputs previously defined and add
# results to the hdf5 file. The visualisation of the output may be done
# with the vview command.
with Hdf5(mode='r+') as io:

    # By default earth gravity is applied and the units are those
    # of the International System of Units.
    # Because of fixed collision margins used in the collision detection,
    # sizes of small objects may need to be expressed in cm or mm.

    # print(pydoc.render_doc(io.run, "Help on %s"))

    io.run(with_timer=False,
           options=options,
           t0=0,
           T=20,
           h=0.005,
           theta=0.50001,
           Newton_max_iter=1,
           set_external_forces=None,
           solver=Numerics.SICONOS_FRICTION_3D_NSGS,
           itermax=10000,
           tolerance=1e-8,
           numerics_verbose=False,
           output_frequency=None)
