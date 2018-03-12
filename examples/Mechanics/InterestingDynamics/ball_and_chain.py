#!/usr/bin/env python

# Various object types sliding, rolling, and sitting still.

from siconos.mechanics.collision.tools import Contactor
from siconos.io.mechanics_io import Hdf5
import siconos.numerics as Numerics
from siconos.mechanics.collision.convexhull import ConvexHull
import numpy as np

# Shape parameters of a single link
num_parts = 8
radius = 0.1
width = 0.03
girth = 0.03
length = 0.2

# Chain and ball parameters
num_links = 5
ball_radius = 0.5

# Creation of the hdf5 file for input/output
with Hdf5() as io:

    # Definition of a half-torus as a composition of convex shapes.
    all_pts = []
    s = np.pi/(num_parts-1)/2
    center = radius-girth
    for i, d in enumerate([-1, 1]):
        for j, r in enumerate(np.linspace(0, np.pi, num_parts)):
            pts = np.array([
                [ width, radius*np.sin(r-s)*d+d*length, radius*np.cos(r-s)*d],
                [-width, radius*np.sin(r-s)*d+d*length, radius*np.cos(r-s)*d],
                [ width, radius*np.sin(r+s)*d+d*length, radius*np.cos(r+s)*d],
                [-width, radius*np.sin(r+s)*d+d*length, radius*np.cos(r+s)*d],
                [ width, center*np.sin(r-s)*d+d*length, center*np.cos(r-s)*d],
                [-width, center*np.sin(r-s)*d+d*length, center*np.cos(r-s)*d],
                [ width, center*np.sin(r+s)*d+d*length, center*np.cos(r+s)*d],
                [-width, center*np.sin(r+s)*d+d*length, center*np.cos(r+s)*d],
            ])
            all_pts += list(pts)
            io.add_convex_shape('Chainlink%02d' % (i*num_parts+j), pts)

    # connector: using the 16 points closest to center, create two
    # hulls for upper and lower part
    all_pts = np.array(all_pts)
    connect = all_pts[np.argsort((all_pts[:, 1]-0)**2)[:16]]
    io.add_convex_shape('Chainlink%02d' % (num_parts*2+0),
                        connect[connect[:, 2] > 0])
    io.add_convex_shape('Chainlink%02d' % (num_parts*2+1),
                        connect[connect[:, 2] < 0])

    # computation of inertia and volume of all points
    ch = ConvexHull(all_pts)
    link_inertia, volume = ch.inertia(ch.centroid())

    # computation of inertia and volume of all points including
    # extrema of the ball
    ball_pos = np.array([0, ball_radius+length*2/3, 0])
    ch = ConvexHull(np.vstack((all_pts,
                               [ball_pos + [0, ball_radius, 0],
                                ball_pos + [0, 0,  ball_radius],
                                ball_pos + [0, 0, -ball_radius],
                                ball_pos + [ball_radius, 0, 0],
                                ball_pos + [-ball_radius, 0, 0]])))
    ball_inertia, volume = ch.inertia(ch.centroid())

    # ball at the end of the chain
    io.add_primitive_shape('Ball', 'Sphere', [ball_radius])

    chainlink = [Contactor('Chainlink%02d' % i) for i in range(num_parts*2+2)]
    ball = [Contactor('Ball', relative_translation=ball_pos)]
    mass = 0.1
    initvel = 0.0
    inertia = link_inertia
    for i in range(num_links):
        io.add_object('link%02d' % (i*2+0), chainlink,
                      translation=[0, 0, 10+(i*2+0)*length*2],
                      orientation=[(1, 0, 0), np.pi/2],
                      velocity=[0, 0, 0, 0, 0, 0],
                      mass=mass*(i > 0), inertia=inertia)
        # Last link has the ball
        if (i+1) == num_links:
            # Ball is 10x heavier, and give slide sideways initial
            # velocity to cause collapse of the chain.
            chainlink += ball
            mass = 1
            initvel = 0.1
            inertia = ball_inertia
        io.add_object('link%02d' % (i*2 + 1), chainlink,
                      translation=[0, 0, 10 + (i*2 + 1)*length*2],
                      orientation=[(1, 1, 1), np.pi*2/3],
                      velocity=[initvel, 0, 0, 0, 0, 0],
                      mass=mass, inertia=inertia)

    # Definition of the ground shape
    io.add_primitive_shape('Ground', 'Box', (10, 10, 1))

    # the ground object made with the ground shape. As the mass is
    # not given, it is a static object only involved in contact
    # detection.
    io.add_object('ground', [Contactor('Ground')],
                  translation=[0, 0, 0])

    # Definition of a non smooth law. As no group ids are specified it
    # is between contactors of group id 0.
    io.add_Newton_impact_friction_nsl('contact', mu=0.03, e=0.0)

# Run the simulation from the inputs previously defined and add
# results to the hdf5 file. The visualisation of the output may be done
# with the vview command.
with Hdf5(mode='r+') as io:

    io.run(with_timer=False,
           t0=0,
           T=20,
           h=0.001,
           theta=0.50001,
           Newton_max_iter=1,
           set_external_forces=None,
           solver=Numerics.SICONOS_FRICTION_3D_NSGS,
           itermax=10000,
           tolerance=1e-8,
           numerics_verbose=False,
           output_frequency=None)
