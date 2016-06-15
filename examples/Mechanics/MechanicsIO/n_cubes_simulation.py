#!/usr/bin/env python

from siconos.io.mechanics_io import Hdf5
import siconos.numerics as Numerics

# Run the simulation from the inputs previously defined and add
# results to the hdf5 file. The visualisation of the output may be done
# with the vview command.


nstep=20000
step=0.0005
with Hdf5(mode='r+',io_filename='n_cubes_scene.hdf5') as io:

    # By default earth gravity is applied and the units are those
    # of the International System of Units.
    # Because of fixed collision margins used in the collision detection,
    # sizes of small objects may need to be expressed in cm or mm.
    io.run(with_timer=False,
           time_stepping=None,
           space_filter=None,
           body_class=None,
           shape_class=None,
           face_class=None,
           edge_class=None,
           gravity_scale=1,
           t0=0,
           T=nstep*step,
           h=step,
           multipoints_iterations=True,
           theta=0.50001,
           Newton_max_iter=1,
           set_external_forces=None,
           solver=Numerics.SICONOS_FRICTION_3D_NSGS,
           itermax=100,
           tolerance=1e-4,
           numerics_verbose=False,
           output_frequency=100)
