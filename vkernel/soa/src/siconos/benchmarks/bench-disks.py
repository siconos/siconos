# Usage: disks-bench <backend> <problem size>

import siconos.numerics as sn
import siconos.simulation as sim

from math import log
import sys

from siconos.mechanics.collision.tools import Contactor
from siconos.mechanics.collision.bullet import SiconosBulletOptions

from siconos.io.mechanics_run import MechanicsHdf5Runner,\
    MechanicsHdf5Runner_run_options,\
    RunnerConfig
import nonos

backend = str(sys.argv[1])
runner_config = RunnerConfig(backend)

N = int(sys.argv[2])

disk_radius = 1

io_filename = 'bench-disks-{}-{}.hdf5'.format(backend,N)

with MechanicsHdf5Runner(io_filename=io_filename, config=runner_config) as io:

    io.add_primitive_shape('DiskR', 'Disk', [disk_radius],
                           insideMargin=1, outsideMargin=0.)

    base = 2*N
    dbase = base
    offset = base * 0.
    vinit = base * 0.
    h1 = 4
    h2 = 7

    l1 = 9
    l2 = 11

    io.add_primitive_shape('Ground-1',
                           'Segment', (-base, h1*base+dbase,
                                       0, h1*base))

    io.add_primitive_shape('Ground-2',
                           'Segment', (0, h1*base,
                                       l1*base, 0))

    io.add_primitive_shape('Ground-3',
                           'Segment', (l1*base, 0,
                                       l2*base, 0))

    io.add_primitive_shape('Wall-1',
                           'Segment', (l2*base, 0, l2*base, h2*base))

    io.add_primitive_shape('Wall-2',
                           'Segment', (-base, h1*base+dbase, -base, h2*base))

    io.add_object('ground-1', [Contactor('Ground-1')], translation=[0,0])
    io.add_object('ground-2', [Contactor('Ground-2')], translation=[0,0])
    io.add_object('ground-3', [Contactor('Ground-3')], translation=[0,0])
    io.add_object('wall-1', [Contactor('Wall-1')], translation=[0,0])
    io.add_object('wall-2', [Contactor('Wall-2')], translation=[0,0])

    for i in range(N):
        for j in range(N):
            io.add_object('disk{}-{}'.format(i,j), [Contactor('DiskR')],
                          translation=[-base+2*i+1, h1*base+2*j+1+dbase + offset],
                          orientation=[0], velocity=[vinit, 0, 0], mass=1, inertia=0.5)

    io.add_Newton_impact_friction_nsl('contact', mu=0.5, e=0)


bullet_options = SiconosBulletOptions()

print("###")
print(dir(bullet_options))
print("###")
#bullet_options.dimension = 2
bullet_options.worldScale = 1.0
bullet_options.perturbationIterations = 1
#bullet_options.contactBreakingThreshold = 1
bullet_options.minimumPointsPerturbationThreshold = 1

options = sn.solver_options_create(sn.solver_ids.SICONOS_FRICTION_2D_NSGS)
options.iparam[sn.params.SICONOS_IPARAM_MAX_ITER] = 100
options.dparam[sn.params.SICONOS_DPARAM_TOL] = 1e-2
options.iparam[sn.params.SICONOS_FRICTION_3D_NSGS_FREEZING_CONTACT] = 10

run_options=MechanicsHdf5Runner_run_options()
run_options['t0']=0
run_options['T']=0.01
run_options['h']=0.005
run_options['gravity_scale'] = 1/N


run_options['bullet_options']=bullet_options
run_options['solver_options']=options

run_options['constraint_activation_threshold']=1e-05
#run_options['Newton_options']=sk.SICONOS_TS_LINEAR

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
run_options['time_stepping']=None


with MechanicsHdf5Runner(io_filename=io_filename, mode='r+',
                         config=runner_config) as io:

        io.run(run_options)

