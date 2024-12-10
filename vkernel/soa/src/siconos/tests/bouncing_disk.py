import siconos.numerics as sn
import siconos.kernel as sk

import sys

from siconos.mechanics.collision.tools import Contactor
from siconos.mechanics.collision.bullet import SiconosBulletOptions
from siconos.io.mechanics_run import MechanicsHdf5Runner, set_backend
from math import pi

backend = str(sys.argv[1])
set_backend(backend)

disk_radius = 1

with MechanicsHdf5Runner() as io:

    io.add_primitive_shape('DiskR', 'Disk', [disk_radius])


#    io.add_primitive_shape('Ground-1',
#                           'Segment', (-10, 0,
#                                       10, 0.))

    io.add_primitive_shape('Ground-2',
                           'Box2d', (20, 1), insideMargin=0., outsideMargin=0.)

#    io.add_object('ground-1', [Contactor('Ground-1')], translation=[0,0], orientation=[pi/4])

    io.add_object('ground-2', [Contactor('Ground-2')], translation=[0,0], orientation=[pi/4])

    io.add_object('disk-1', [Contactor('DiskR')],
                  translation=[0, 10],
                  orientation=[0], velocity=[0, 0, 0], mass=1, inertia=0.5)

    io.add_Newton_impact_friction_nsl('contact', mu=0.5, e=0)

bullet_options = SiconosBulletOptions()
bullet_options.worldScale = 1.0
bullet_options.perturbationIterations = 1
bullet_options.minimumPointsPerturbationThreshold = 1

options = sk.solver_options_create(sn.SICONOS_FRICTION_2D_NSGS)
options.iparam[sn.SICONOS_IPARAM_MAX_ITER] = 20
options.dparam[sn.SICONOS_DPARAM_TOL] = 1e-4

def noforces(body):
    pass

with MechanicsHdf5Runner(mode='r+') as io:

        io.run(with_timer=True,
               bullet_options=bullet_options,
               gravity_scale=1,
               t0=0,
               T=10,
               h=0.005,
               theta=0.50001,
               Newton_max_iter=1,
               set_external_forces=None,
               solver_options=options,
               numerics_verbose=True,
               numerics_verbose_level=1,
               output_contact_forces=True,
               output_frequency=None)

