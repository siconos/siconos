#!/usr/bin/env python

#
# A circle with two disks inside under earth gravity
#

from siconos.mechanics.collision.tools import Contactor
from nonos.mechanics_run import MechanicsHdf5Runner

import nonos
import siconos.numerics as sn
import siconos.kernel as sk

import siconos

import numpy
from math import sqrt

nonos.mechanics_run.set_backend('vnative')

disk_radius = 2
circle_radius = 10


def make_input():

    # Creation of the hdf5 file for input/output
    with MechanicsHdf5Runner() as io:

        # Definition of a disk radius 1
        io.add_primitive_shape('DiskR', 'Disk', [disk_radius])

        io.add_primitive_shape('CircleR', 'Circle', [circle_radius])

        # Definition of the ground shape
        io.add_primitive_shape('Ground', 'Line', (0, 30, 0))

    #    io.add_primitive_shape('Wall1', 'Line', (1, 0, -20))

    #    io.add_primitive_shape('Wall2', 'Line', (1, 0, 20))

        # Definition of a non smooth law. As no group ids are specified it
        # is between contactors of group id 0.
        io.add_Newton_impact_friction_nsl('contact', mu=0.3, e=0)

        # The disk object made with an unique Contactor : the Disk shape.
        # As a mass is given, it is a dynamic system involved in contact
        # detection and in the simulation.  With no group id specified the
        # Contactor belongs to group 0
        io.add_object('disk0', [Contactor('DiskR')],
                      translation=[-(circle_radius-disk_radius), circle_radius],
                      orientation=[0], velocity=[0, 0, 0], mass=10)

        io.add_object('disk1', [Contactor('DiskR')],
                      translation=[(circle_radius-disk_radius), circle_radius],
                      orientation=[0], velocity=[0, 0, 0], mass=10)

        io.add_object('circle', [Contactor('CircleR')],
                      translation=(0, circle_radius),
                      orientation=[0], velocity=[0, 0, 0], mass=1)

        # the ground object made with the ground shape. As the mass is
        # not given, it is a static object only involved in contact
        # detection.
        io.add_object('ground', [Contactor('Ground')],
                      translation=[0, 0])

    #    io.add_object('wall1', [Contactor('Wall1')],
    #                  translation=[0, 0])

    #    io.add_object('wall2', [Contactor('Wall2')],
    #                  translation=[0, 0])


# Run the simulation from the inputs previously defined and add
# results to the hdf5 file. The visualisation of the output may be done
# with the vview command.

# LEMKE failure when mu=0
options = sk.solver_options_create(sn.SICONOS_FRICTION_2D_NSGS)
options.iparam[sn.SICONOS_IPARAM_MAX_ITER] = 100000
options.dparam[sn.SICONOS_DPARAM_TOL] = 1e-12


def run():

    with MechanicsHdf5Runner(mode='r+') as io:

        # By default earth gravity is applied and the units are those
        # of the International System of Units.
        # Because of fixed collision margins used in the collision detection,
        # sizes of small objects may need to be expressed in cm or mm.
        io.run(with_timer=False,
               gravity_scale=1,
               t0=0,
               T=2,
               h=0.005,
               theta=0.50001,
               Newton_max_iter=1000,
               set_external_forces=None,
               solver_options=options,
               numerics_verbose=False,
               output_contact_forces=False,
               output_frequency=None)


def check():

    with MechanicsHdf5Runner(
            interaction_manager=nonos.disks.interaction_manager,
            mode='r') as io:

        positions = io.dynamic_data()
        velocities = io.velocities_data()

        times, indices = numpy.unique(positions[:, 0], return_index=True)
        t = 1.95
        id_t = max(0, numpy.searchsorted(times, t, side='right') - 1)

        if id_t < len(indices)-1:
            id_t_m = list(range(indices[id_t],
                                indices[id_t+1]))
        else:
            id_t_m = [indices[id_t]]

        disk1_vidx = numpy.argwhere(velocities[id_t_m, 1] == 1)
        disk2_vidx = numpy.argwhere(velocities[id_t_m, 1] == 2)
        circle_vidx = numpy.argwhere(velocities[:, 1] == 3)

        circle_vidx= circle_vidx.reshape(circle_vidx.shape[0])

        disk1_evx = velocities[id_t_m, 2][disk1_vidx]
        disk1_evy = velocities[id_t_m, 3][disk1_vidx]

        disk2_evx = velocities[id_t_m, 2][disk2_vidx]
        disk2_evy = velocities[id_t_m, 3][disk2_vidx]

        # disks velocities at the end are 0
        assert disk2_evx**2 + disk2_evy**2 < 1e-10
        assert disk1_evx**2 + disk1_evy**2 < 1e-10

        circle_vx = velocities[circle_vidx, 2]
        circle_vy = velocities[circle_vidx, 3]

        # circle velocities should always be near 0
        print('max circle velocity:', max(circle_vx**2 + circle_vy**2))

        assert max(circle_vx**2 + circle_vy**2) < 1e-10

        disk1_idx = numpy.argwhere(positions[:, 1] == 1)
        disk2_idx = numpy.argwhere(positions[:, 1] == 2)
        circle_idx = numpy.argwhere(positions[:, 1] == 3)

        disk1_idx= disk1_idx.reshape(disk1_idx.shape[0])
        disk2_idx= disk2_idx.reshape(disk2_idx.shape[0])
        circle_idx= circle_idx.reshape(circle_idx.shape[0])

        disk1_x = positions[disk1_idx, 2]
        disk1_y = positions[disk1_idx, 3]

        disk2_x = positions[disk2_idx, 2]
        disk2_y = positions[disk2_idx, 3]

        circle_x = positions[circle_idx, 2]
        circle_y = positions[circle_idx, 3]

        print('max disk1 dist:', sqrt(max((disk1_x-circle_x)**2 +
                                          (disk1_y-circle_y)**2))
              - (circle_radius - disk_radius))
        print('max disk2 dist:', sqrt(max((disk2_x-circle_x)**2 +
                                          (disk2_y-circle_y)**2))
              - (circle_radius - disk_radius))

        assert ((sqrt(max((disk1_x-circle_x)**2 + (disk1_y-circle_y)**2)) -
                 (circle_radius - disk_radius)) < 1e-4)
        assert ((sqrt(max((disk2_x-circle_x)**2 + (disk2_y-circle_y)**2)) -
                 (circle_radius - disk_radius)**2) < 1e-4)
        assert sqrt(max(circle_x**2 + (circle_y-circle_radius)**2)) < 1e-10


def test_native_collision():
    make_input()
    run()
    check()

if __name__ == "__main__":
    test_native_collision()
