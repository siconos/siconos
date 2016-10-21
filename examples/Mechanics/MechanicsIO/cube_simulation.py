
from siconos.io.mechanics_io import Hdf5

import siconos.io.mechanics_io
siconos.io.mechanics_io.set_implementation('original')
siconos.io.mechanics_io.set_backend('bullet')

# Run the simulation from the inputs previously defined and add
# results to the hdf5 file. The visualisation of the output may be done
# with the vview command.
with Hdf5(mode='r+',io_filename='cube_scene.hdf5') as io:

    # By default earth gravity is applied and the units are those
    # of the International System of Units.
    # Because of fixed collision margins used in the collision detection,
    # sizes of small objects may need to be expressed in cm or mm.
    io.run(output_frequency=100)
