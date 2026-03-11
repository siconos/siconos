import siconos.numerics as sn
import siconos.simulation as sim

from siconos.mechanics.collision.tools import Contactor
from siconos.mechanics.collision.bullet import SiconosBulletOptions

from siconos.io.mechanics_run import MechanicsHdf5Runner,\
    MechanicsHdf5Runner_run_options,\
    RunnerConfig

import nonos

import sys
import numpy
from math import sqrt

backend = str(sys.argv[1])
runner_config = RunnerConfig(backend)

io_filename = f'one-disk-{backend}.hdf5'

disk_radius = 1

with MechanicsHdf5Runner(io_filename=io_filename, config=runner_config) as io:

    io.add_primitive_shape('DiskR', 'Disk', [disk_radius])

    io.add_primitive_shape('Ground', 'Segment', (-30, 1, 30, 1))

    io.add_Newton_impact_friction_nsl('contact', mu=0.3, e=0)

    io.add_object('disk0', [Contactor('DiskR')],
                      translation=[0, 3*disk_radius],
                      orientation=[0], velocity=[0, 0, 0], mass=1, inertia=0.5)

    io.add_object('ground', [Contactor('Ground')],
                  translation=[0, 0])

bullet_options = SiconosBulletOptions()
bullet_options.worldScale = 1.0
bullet_options.perturbationIterations = 1
#bullet_options.contactBreakingThreshold = 1
bullet_options.minimumPointsPerturbationThreshold = 1

options = sn.solver_options_create(sn.solver_ids.SICONOS_FRICTION_2D_NSGS)
options.iparam[sn.params.SICONOS_IPARAM_MAX_ITER] = 100
options.dparam[sn.params.SICONOS_DPARAM_TOL] = 1e-2
options.iparam[sn.params.SICONOS_NSGS_FREEZING_CONTACT] = 10

run_options=MechanicsHdf5Runner_run_options()
run_options['t0']=0
run_options['T']=1
run_options['h']=0.05

run_options['bullet_options']=bullet_options
run_options['solver_options']=options

run_options['constraint_activation_threshold']=1e-05
run_options['Newton_options']=sim.LINEAR

run_options['skip_last_update_output']=True
run_options['skip_reset_lambdas']=True
#run_options['osns_assembly_type']= sk.REDUCED_DIRECT

#run_options['osns_assembly_type']= sk.GLOBAL_REDUCED
#run_options['osi']= sk.MoreauJeanGOSI
#run_options['skip_last_update_input']=True
#if performance_verbose:
#    run_options['verbose']=True
#    run_options['with_timer']=True
#    run_options['explode_Newton_solve']=True
#    run_options['explode_computeOneStep']=True

run_options['violation_verbose']=True
run_options['output_frequency']=1
run_options['with_timer']=True
#run_options['time_stepping']=None
#run_options['output_contact_info']=True
#run_options['output_contact_work']=False # failure for bullet/features-eigen2

#options = sn.solver_options_create(sn.solver_ids.SICONOS_FRICTION_2D_LEMKE)



with MechanicsHdf5Runner(io_filename=io_filename, mode='r+',
                         config=runner_config) as io:

        io.run(run_options)
