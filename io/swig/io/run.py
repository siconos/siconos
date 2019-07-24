#!/usr/bin/env @PYTHON_EXECUTABLE@
"""
Description: Run a pre-generated Siconos mechanics-IO HDF5 simulation file.
"""

# Lighter imports before command line parsing
from __future__ import print_function
import argparse

parser = argparse.ArgumentParser(
    description = __doc__,
    epilog = """This script only provides a basic interface for the most common
    simulation options.  For more complex options, or to specify
    behaviour such as controllers, you must create a custom simulation
    script.  If a partially-completed simulation is found, it will be
    continued from the last time step until T.

    Note that most example scripts do not use this program, and simply
    define and then run the simulation in the same script.""")

parser.add_argument('file', metavar='filename', type=str, nargs=1,
                    help = 'simulation file (HDF5)')
parser.add_argument('-T', metavar='time', type=float,
                    help = 'time in seconds to run until (default T=1 second)',
                    default=1)
parser.add_argument('-p', '--period', metavar='period', type=float,
                    help = 'time in seconds between frames (default p=5e-3)',
                    default=5e-3)
parser.add_argument('-e','--every', metavar='interval', type=int,
                    help = 'output every nth frame (default=1)', default=1)
parser.add_argument('-f','--frequency', metavar='Hz', type=float,
                    help = 'alternative to -p, specify simulation frequency in Hz'+
                    ' (default=200 Hz)')
parser.add_argument('-V','--version', action='version',
                    version='@SICONOS_VERSION@')

args = parser.parse_args()

if args.frequency is not None:
    args.p = 1.0 / args.frequency

# Heavier imports after command line parsing
from siconos.io.mechanics_run import MechanicsHdf5Runner

# Run the simulation from the inputs previously defined and add
# results to the hdf5 file. The visualisation of the output may be done
# with the vview command.
with MechanicsHdf5Runner(mode='r+',io_filename=args.file[0]) as io:

    # By default earth gravity is applied and the units are those
    # of the International System of Units.
    # Because of fixed collision margins used in the collision detection,
    # sizes of small objects may need to be expressed in cm or mm.
    io.run(output_frequency=args.every, T=args.T, h=args.period)
