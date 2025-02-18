from siconos.mechanics.collision.tools import Contactor
from siconos.io.mechanics_run import MechanicsHdf5Runner, \
    MechanicsHdf5Runner_run_options, set_backend

import nonos
import siconos.numerics as sn
import siconos.kernel as sk

import siconos

import numpy
from math import sqrt
from random import random

set_backend('vnative')

disk_radius = 1
N = 10

with MechanicsHdf5Runner() as io:

    io.add_primitive_shape('DiskR1', 'Disk', [disk_radius])
    io.add_primitive_shape('DiskR2', 'Disk', [0.1*disk_radius])

    io.add_Newton_impact_friction_nsl('contact', mu=0.3, e=0)

#    io.add_primitive_shape('Ground', 'Segment', [-10,-10,10,-10])
    for i in range(N):
        for j in range(N):
            io.add_object('disk-{}-{}'.format(i,j), [Contactor('DiskR2')],
                          translation=[0.2*i, 3*disk_radius+0.2*j],
                          orientation=[0], velocity=[0, 0, 0], mass=1, inertia=0.5)

    io.add_object('fixed-disk1', [Contactor('DiskR1')], translation=[0, 0])

    io.add_object('fixed-disk2', [Contactor('DiskR1')], translation=[2, -5])

#    io.add_object('ground' , [Contactor('Ground')], translation=[0, 0])

options = sk.solver_options_create(sn.SICONOS_FRICTION_2D_NSGS)
options.iparam[sn.SICONOS_IPARAM_MAX_ITER] = 10
options.dparam[sn.SICONOS_DPARAM_TOL] = 1e-3

class recirculation_start_hook():

    def __init__(self):
        pass

    def initialize(self, io):
        self._io = io

    def call(self, step):
        positions  = self._io._io.positions(self._io._nsds)
        y = positions[:,2]
        maxy=numpy.max(y)
        ds_idy = numpy.nonzero(y < -6)

        for i,k in enumerate(ds_idy[0]):

            nds = int(positions[k, 0])
            ds = self._io._nsds.dynamicalSystem(nds)
#            print('ds {} is under 0! : {}'.format(nds, ds.handle().q()))

            jj = nds

            ds.setQ0Ptr([(random()-0.5)/100, maxy + 0.2*(i+1), 0.])
            ds.setVelocity0Ptr([0., 0., 0.])
            ds.resetToInitialState()
            ds.swapInMemory()


run_options=MechanicsHdf5Runner_run_options()

run_options['t0']=0
run_options['T']=10
run_options['h']=0.005
run_options['theta']=0.50001
run_options['solver_options']=options
run_options['start_run_iteration_hook']=recirculation_start_hook()

with MechanicsHdf5Runner(mode='r+') as io:
        io.run(run_options)
