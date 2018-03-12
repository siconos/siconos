#!/usr/bin/env python

from siconos.io.mechanics_io import Hdf5
import siconos.numerics as Numerics
import siconos.kernel as Kernel

import chute
import rocas
import random

random.seed(0)

unscaled_cube_size = 0.1
unscaled_plan_thickness = unscaled_cube_size
unscaled_density = 2500

scale = 1.0/unscaled_cube_size*10.0

cube_size = unscaled_cube_size*scale
plane_thickness = unscaled_plan_thickness*scale
density =  unscaled_density/(scale**3)

box_height = 3.683
box_length = 6.900
box_width  = 3.430

plane_thickness = 0.2

with Hdf5(mode='w') as io:
    ch = chute.create_chute(io, box_height = box_height,
                            box_length = box_length,
                            box_width = box_width,
                            plane_thickness = plane_thickness,
                            scale = 1, trans = [-0.6, -1.8, -1])

    rcs = rocas.create_rocas(io, n_layer=200, n_row=2, n_col=16,
                             x_shift=2.0, roca_size=0.1, top=3,
                             rate=0.02, density=density)

    io.add_Newton_impact_friction_nsl('contact', mu=1.0, e=0.01)

step=20000
hstep=1e-4
itermax=10000
dump_probability = .02
theta=1.0
tolerance=1e-03
import os 
from os import path
if not os.path.exists('Chute'):
    os.mkdir('Chute')


solver=Numerics.SICONOS_GLOBAL_FRICTION_3D_ADMM
fileName = "./Chute/Chute"
title = "Chute"
description = """
Chute with 6400 polyhedra with Bullet collision detection
Moreau TimeStepping: h={0}, theta = {1}
One Step non smooth problem: {2}, maxiter={3}, tol={4}
""".format(hstep, theta, Numerics.solver_options_id_to_name(solver),
           itermax,
           tolerance)
mathInfo = ""

from siconos.io.FrictionContactTrace import FrictionContactTraceParams
friction_contact_trace_params = FrictionContactTraceParams(dump_itermax=9000, dump_probability=None,
                                                           fileName=fileName, title =title,
                                                           description = description, mathInfo= mathInfo)

    
with Hdf5(mode='r+', collision_margin=0.01) as io:
    # By default earth gravity is applied and the units are those
    # of the International System of Units.
    # Because of fixed collision margins used in the collision detection,
    # sizes of small objects may need to be expressed in cm or mm.
    io.run(gravity_scale=1.0/scale,
           t0=0,
           T=step*hstep,
           h=hstep,
           multipoints_iterations=True,
           theta=1.0,
           Newton_max_iter=1,
           solver=Numerics.SICONOS_FRICTION_3D_NSGS,
           itermax=10000,
           tolerance=1e-3,
           output_frequency=10,
           osi=Kernel.MoreauJeanGOSI,
           friction_contact_trace=True,
           friction_contact_trace_params=friction_contact_trace_params)
