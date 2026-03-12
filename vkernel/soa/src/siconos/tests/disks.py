import siconos.numerics as sn
import siconos.kernel as sk

import numpy
from math import sqrt

import random


from siconos.mechanics.collision.tools import Contactor
try:
    from nonos.mechanics_run import MechanicsHdf5Runner
    import nonos
    nonos.mechanics_run.set_backend('vnative')
except Exception:
    from siconos.io.mechanics_run import MechanicsHdf5Runner, \
    MechanicsHdf5Runner_run_options, set_backend
    sn.numerics_set_verbose(True)
    set_backend('native')


disk_radius = 1

with MechanicsHdf5Runner() as io:

    io.add_primitive_shape('DiskR', 'Disk', [disk_radius])

    io.add_primitive_shape('Ground-0',
                           'Segment', (-10, 10, 0, 0))

    io.add_primitive_shape('Ground-1',
                           'Segment', (0, 0, 10, 10))

    io.add_Newton_impact_friction_nsl('contact', mu=0.5, e=0)

    # io.add_object('disk0', [Contactor('DiskR')],
    #               translation=[8, 8],
    #               orientation=[0], velocity=[0, 0, -5], mass=1)

    # io.add_object('disk1', [Contactor('DiskR')],
    #               translation=[3, 5],
    #               orientation=[0], velocity=[0, 0, 0], mass=1)

    for i in range(2):
        io.add_object('ground-{}'.format(i), [Contactor('Ground-{}'.format(i))],
                      translation=[0, 0])

    # io.add_object('ground2', [Contactor('Ground2')],
    #               translation=[0, 0])

    # io.add_object('ground3', [Contactor('Ground3')],
    #               translation=[0, 0])

    io.add_object('disk{}-{}'.format(0,0), [Contactor('DiskR')],
                  translation=[0,5], orientation=[0], velocity=[0,0,0], mass=1)


    io.add_object('disk{}-{}'.format(0,1), [Contactor('DiskR')],
                  translation=[0,2.01], orientation=[0], velocity=[0,0,0], mass=1)

    # N = 1
    # for i in range(N):
    #     for j in range(N):
    #         io.add_object('disk{}-{}'.format(i,j), [Contactor('DiskR')],
    #                       translation=[2*i-N, 2*j+2+.333333334*N],
    #                       orientation=[0], velocity=[0, 0, 0], mass=1)


    # io.add_object('disk2', [Contactor('DiskR')],
    #               translation=[0.5, 4+3*disk_radius],
    #               orientation=[0], velocity=[0, 0, 0], mass=1)


options = sk.solver_options_create(sn.SICONOS_FRICTION_2D_NSGS)
options.iparam[sn.SICONOS_IPARAM_MAX_ITER] = 20
options.dparam[sn.SICONOS_DPARAM_TOL] = 1e-4

def noforces(body):
    pass

with MechanicsHdf5Runner(mode='r+') as io:

        io.run(with_timer=False,
               gravity_scale=1,
               t0=0,
               T=20,
               h=0.005,
               theta=0.50001,
               Newton_max_iter=1000,
               set_external_forces=None,
               solver_options=options,
               numerics_verbose=True,
               numerics_verbose_level=1,
               output_contact_forces=True,
               output_frequency=None)
