#!/usr/bin/env python

from siconos.io.mechanics_io import Hdf5
import siconos.numerics as Numerics

#import chute_con_vibrator_rear_up as chute_con_vibradores
#import chute_con_vibrador_bottom as chute_con_vibradores
#import chute_con_vibrador_middle as chute_con_vibradores
import chute_con_vibrador_top as chute_con_vibradores

import rocas
import random
import sys
import numpy

if (len(sys.argv) < 2):
    dist = 'uniform'
    mu = 0.1
else:
    dist = sys.argv[1]
    mu = sys.argv[2]




if not dist in ['uniform', 'double', 'exp']:
    print("dist = [uniform | double | exp]")
    sys.exit(1)
if float(mu) < 0.1 or float(mu) > 1.5:
    print("mu = [0.1 .. 1.5]")
    sys.exit(1)

fn = 'chute_con_rocas_y_vibradores-{0}-mu-{1}.hdf5'.format(dist,mu)

random.seed(0)
numpy.random.seed(0)

cube_size = 0.1
plan_thickness = cube_size
density = 2500

box_height = 3.683
box_length = 6.900
box_width  = 3.430

plane_thickness = 0.2

with Hdf5(mode='w', io_filename=fn) as io:

    io.addNewtonImpactFrictionNSL('contact_rocas', mu=float(mu), e=0.01,collision_group1=0, collision_group2=0)
    io.addNewtonImpactFrictionNSL('contact_chute_rocas', mu=float(mu), e=0.01,collision_group1=1, collision_group2=0)
    io.addNewtonImpactFrictionNSL('contact_vibratores_rocas', mu=float(mu), e=0.01,collision_group1=2, collision_group2=0)


    
    ch = chute_con_vibradores.create_chute(io, box_height = box_height,
                                           box_length = box_length,
                                           box_width = box_width,
                                           plane_thickness = plane_thickness,
                                           scale = 1, trans = [-0.6, -1.8, -1])
    
    rcs = rocas.create_rocas(io, n_layer=20, n_row=2, n_col=5,
                             x_shift=1.5, roca_size=0.13, top=3,
                             rate=0.2, density=density,
                             distribution = (dist, 0.12))

 
    
    

hstep =1e-3
with Hdf5(mode='r+', io_filename=fn) as io:
    io.run(with_timer=False,
           t0=0,
           T=45, #100*hstep,
           h=hstep,
           multipoints_iterations=True,
           theta=0.50001,
           Newton_max_iter=1,
           set_external_forces=None,
           solver=Numerics.SICONOS_FRICTION_3D_NSGS,
           itermax=1000,
           tolerance=1e-4,
           numerics_verbose=False,
           output_frequency=10)
