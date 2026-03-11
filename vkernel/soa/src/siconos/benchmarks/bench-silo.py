import siconos.numerics as sn
import siconos.simulation as sim
import sys
import os
from siconos.mechanics.collision.tools import Contactor
from siconos.mechanics.collision.bullet import SiconosBulletOptions
from siconos.io.mechanics_run import MechanicsHdf5Runner, MechanicsHdf5Runner_run_options, RunnerConfig
from siconos.io.FrictionContactTrace import FrictionContactTraceParams
import random
import numpy as np
from math import sqrt

backend = str(sys.argv[1])
runner_config = RunnerConfig(backend)

sqrtN = int(sys.argv[2])
itermax = int(sys.argv[3])
tolerance = float(sys.argv[4])
hstep = float(sys.argv[5])
neighborhood_radius_factor = float(sys.argv[6])
neighborhood_min_radius_factor = float(sys.argv[7])


N = sqrtN*sqrtN

disk_radius = 10/sqrtN
io_filename = f'bench-silo-disks-{backend}-{sqrtN}-{itermax}-{tolerance}-{hstep}-{neighborhood_radius_factor}-{neighborhood_min_radius_factor}.hdf5'


sqrtNp = 10

silo_offset = sqrtNp

silo_width = 2 * sqrtNp
silo_height = 2 * sqrtNp + silo_offset + disk_radius

box_width = silo_width * 3
box_height = silo_height // 2

disk_mass = 1
disk_inertia = 1

if not os.path.exists(io_filename):
    with MechanicsHdf5Runner(io_filename=io_filename, config=runner_config) as io:


        radii = disk_radius * np.array([1., 0.95, 0.9, 0.85])
        io.add_primitive_shape('DiskR1', 'Disk', [radii[0]], insideMargin=1, outsideMargin=0.)
        io.add_primitive_shape('DiskR2', 'Disk', [radii[1]], insideMargin=1, outsideMargin=0.)
        io.add_primitive_shape('DiskR3', 'Disk', [radii[2]], insideMargin=1, outsideMargin=0.)
        io.add_primitive_shape('DiskR4', 'Disk', [radii[3]], insideMargin=1, outsideMargin=0.)

        disk_contactors = list(zip(radii, [Contactor('DiskR1'),
                                           Contactor('DiskR2'),
                                           Contactor('DiskR3'),
                                           Contactor('DiskR4')]))

        # Silo geometry
        # Left silo wall
        io.add_primitive_shape(
            'SiloWallLeft',
            'Segment',
            (0, silo_offset, 0, silo_height)
        )
        # Right silo wall
        io.add_primitive_shape(
            'SiloWallRight',
            'Segment',
            (silo_width, silo_offset, silo_width, silo_height)
        )

        # Silo bottom wall with a gap in the middle (opening)
        opening_size = silo_width / 4
        # left part of bottom
        io.add_primitive_shape(
            'SiloBottomLeft',
            'Segment',
            (0, silo_offset, (silo_width - opening_size) / 2, 0)
        )
        # right part of bottom
        io.add_primitive_shape(
            'SiloBottomRight',
            'Segment',
            ((silo_width + opening_size) / 2, 0, silo_width, silo_offset)
        )

        # Under the silo, a box to collect disks
        # Box left wall
        io.add_primitive_shape(
            'BoxWallLeft',
            'Segment',
            ((silo_width - box_width) // 2, -box_height, (silo_width - box_width) // 2, 0)
        )
        # Box right wall
        io.add_primitive_shape(
            'BoxWallRight',
            'Segment',
            ((silo_width + box_width) // 2, -box_height, (silo_width + box_width) // 2, 0)
        )
        # Box bottom
        io.add_primitive_shape(
            'BoxBottom',
            'Segment',
            ((silo_width - box_width) // 2, -box_height, (silo_width + box_width) // 2, -box_height)
        )

        # Add static objects
        io.add_object('silo-wall-left', [Contactor('SiloWallLeft')], translation=[0, 0])
        io.add_object('silo-wall-right', [Contactor('SiloWallRight')], translation=[0, 0])
        io.add_object('silo-bottom-left', [Contactor('SiloBottomLeft')], translation=[0, 0])
        io.add_object('silo-bottom-right', [Contactor('SiloBottomRight')], translation=[0, 0])

        io.add_object('box-wall-left', [Contactor('BoxWallLeft')], translation=[0, 0])
        io.add_object('box-wall-right', [Contactor('BoxWallRight')], translation=[0, 0])
        io.add_object('box-bottom', [Contactor('BoxBottom')], translation=[0, 0])

        # Generate positions for disks
        y_offset = silo_height
        disks_per_row = sqrtN
        num_rows = sqrtN

        x_vals = disk_radius + np.arange(disks_per_row) * (2 * disk_radius)
        y_vals = y_offset + disk_radius - np.arange(num_rows) * (2 * disk_radius)

        xx, yy = np.meshgrid(x_vals, y_vals, indexing='xy')

        # Flatten and apply silo width constraint
        xx_flat = xx.ravel()
        yy_flat = yy.ravel()
        mask = xx_flat + disk_radius <= silo_width

        # Stack coordinates and take first N positions
        positions = np.column_stack([xx_flat[mask], yy_flat[mask]])[:N]

        if len(positions) == 0:
            raise ValueError("No valid positions generated! Check silo_width vs disk_radius")

        # Shuffle in-place (numpy version)
        np.random.shuffle(positions)

        # Distribute among different disk types and create aggregates
        n_types = len(disk_contactors)
        base_count = N // n_types
        remainder = N % n_types

        idx = 0
        for i, (radius, contactor) in enumerate(disk_contactors):
            count = base_count + (1 if i < remainder else 0)
            if count == 0:
                continue

            # Slice the pre-shuffled array
            type_positions = positions[idx:idx+count]
            idx += count

            # type_positions is already a float64 numpy array
            translations = type_positions
            orientations = np.zeros(count, dtype=np.float64)
            velocities = np.zeros((count, 3), dtype=np.float64)

            io.add_objects(
                f'granular_disks_{i}',
                [contactor],
                translations=translations,
                orientations=orientations,
                velocities=velocities,
                mass=disk_mass * radius * radius,
                inertia=disk_inertia
            )

        io.add_Newton_impact_friction_nsl('contact', mu=0.5, e=0)

bullet_options = SiconosBulletOptions()
bullet_options.worldScale = 1.0
bullet_options.perturbationIterations = 1
bullet_options.minimumPointsPerturbationThreshold = 1

theta = 0.50001

solver = sn.solver_ids.SICONOS_FRICTION_2D_NSGS

options = sn.solver_options_create(solver)
options.iparam[sn.params.SICONOS_IPARAM_MAX_ITER] = itermax
options.dparam[sn.params.SICONOS_DPARAM_TOL] = tolerance
options.iparam[sn.params.SICONOS_NSGS_FREEZING_CONTACT] = 10

fileName = f"./Silo-{N}-2D-disks"
title = "Silo simulation"
description = f"""
2D disks falling in a box through a silo.
Moreau TimeStepping: h={hstep}, theta = {theta}
One Step non smooth problem: SICONOS_FRICTION_2D_NSGS,
maxiter={itermax}, tol={tolerance},
neighborhood_radius_factor={neighborhood_radius_factor},
neighborhood_min_radius_factor={neighborhood_min_radius_factor}
"""

mathInfo = ""

friction_contact_trace_params = FrictionContactTraceParams(
     dump_itermax=10, dump_probability=None,
     fileName=fileName, title=title,
     description=description, mathInfo=mathInfo)

run_options = MechanicsHdf5Runner_run_options()
run_options['t0'] = 0
run_options['T'] = 15
run_options['h'] = hstep
run_options['gravity_scale'] = 1
run_options['bullet_options'] = bullet_options
run_options['solver_options'] = options
run_options['constraint_activation_threshold'] = 1e-05
#run_options['Newton_options'] = sim.SICONOS_TS_LINEAR
run_options['skip_last_update_output'] = True
run_options['skip_reset_lambdas'] = True
#run_options['osns_assembly_type'] = sim.REDUCED_DIRECT
run_options['violation_verbose'] = True
run_options['output_frequency'] = 0.04 / hstep # 25 fps
run_options['with_timer'] = True
run_options['time_stepping'] = None
#run_options['friction_contact_trace'] = True
#run_options['friction_contact_trace_params'] = friction_contact_trace_params
run_options['output_backup'] = True
run_options['output_backup_frequency'] = 100

import nonos
vnative_options=nonos.bridge.SpaceFilterOptions()
vnative_options.neighborhood_radius = neighborhood_radius_factor * disk_radius
vnative_options.min_radius = disk_radius * neighborhood_min_radius_factor

run_options['vnative_options'] = vnative_options

with MechanicsHdf5Runner(io_filename=io_filename, mode='r+',
                         config=runner_config) as io:
    io.run(run_options)
