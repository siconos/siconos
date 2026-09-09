#import petsc4py

import siconos.numerics as sn
import siconos.simulation as sim

import sys

from siconos.mechanics.collision.tools import Contactor
from siconos.mechanics.collision.bullet import SiconosBulletOptions

from siconos.io.mechanics_run import MechanicsHdf5Runner,\
    MechanicsHdf5Runner_run_options,\
    RunnerConfig
import nonos

shape_filename=sys.argv[1]

backend = 'vnative' # not implemented for bullet
runner_config = RunnerConfig(backend)

disk_radius = .1

with MechanicsHdf5Runner(config=runner_config) as io:

    io.add_primitive_shape('DiskR', 'Disk', [disk_radius])

    io.add_object('disk-1', [Contactor('DiskR')],
                  translation=[0.8, 1.5],
                  orientation=[0], velocity=[0, 0, 0], mass=1, inertia=1)
    
    io.add_primitive_shape('Ground-1',
                           'Segment', (-10, 0,
                                       10,  0))
    
    io.add_shape_data_from_file('Square', shape_filename)
    io.add_object('square', [Contactor('Square')],
                  # density : 2500
                  # young : 1.e11
                  # poisson : 0.3
                  material=(2500,1e11,0.3))

    io.add_Newton_impact_friction_nsl('contact', mu=0.5, e=0)

options = sn.SolverOptions(sn.solver_ids.SICONOS_FRICTION_2D_NSGS)
options.iparam[sn.params.SICONOS_IPARAM_MAX_ITER] = 100
options.dparam[sn.params.SICONOS_DPARAM_TOL] = 1e-2
options.iparam[sn.params.SICONOS_NSGS_FREEZING_CONTACT] = 10
    
with MechanicsHdf5Runner(mode='r+', config=runner_config) as io:

        io.run(with_timer=True,
               t0=0,
               T=2,
               h=0.005,
               theta=0.5001,
               Newton_max_iter=1,
               set_external_forces=None,
               solver_options=options,
               numerics_verbose=True,
               numerics_verbose_level=1,
               output_contact_forces=True,
               output_frequency=None)



        
