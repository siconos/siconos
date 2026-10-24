#!/usr/bin/env @Python_EXECUTABLE@
"""
Description: Run a pre-generated Siconos mechanics-IO HDF5 simulation file.
"""

# Lighter imports before command line parsing
import argparse

parser = argparse.ArgumentParser(
    description=__doc__,
    epilog="""This script only provides a basic interface for the most common
    simulation options.  For more complex options, or to specify
    behaviour such as controllers, you must create a custom simulation
    script.  If a partially-completed simulation is found, it will be
    continued from the last time step until T.

    Note that most example scripts do not use this program, and simply
    define and then run the simulation in the same script.""",
)

parser.add_argument(
    "file", metavar="filename", type=str, nargs=1, help="simulation file (HDF5)"
)
parser.add_argument(
    "-T",
    metavar="time",
    type=float,
    help="time in seconds to run until (default T=1 second)",
    default=1,
)
parser.add_argument(
    "-p",
    "--period",
    metavar="period",
    type=float,
    help="time in seconds between frames (default p=5e-3)",
    default=5e-3,
)
parser.add_argument(
    "-e",
    "--every",
    metavar="interval",
    type=int,
    help="output every nth frame (default=1)",
    default=1,
)
parser.add_argument(
    "-f",
    "--frequency",
    metavar="Hz",
    type=float,
    help="alternative to -p, specify simulation frequency in Hz" + " (default=200 Hz)",
)
parser.add_argument(
    "-b",
    "--backend",
    type=str,
    default="bullet",
    help="collision backend: bullet, native or vnative (vectorized kernel); "
    "default bullet",
)
parser.add_argument(
    "--neighborhood-radius",
    type=float,
    help="vnative: broadphase radius, larger than the largest body extent "
    "plus the contact gap (default: the vnative default)",
)
parser.add_argument(
    "--min-radius",
    type=float,
    help="vnative: smallest body size, sets the sampling of the static "
    "segments (default: the vnative default)",
)
parser.add_argument("-V", "--version", action="version", version="@SICONOS_VERSION@")

args = parser.parse_args()

if args.frequency is not None:
    args.period = 1.0 / args.frequency

# Heavier imports after command line parsing
from siconos.io.mechanics_run import (
    MechanicsHdf5Runner,
    MechanicsHdf5Runner_run_options,
    RunnerConfig,
)

run_options = MechanicsHdf5Runner_run_options()
run_options["output_frequency"] = args.every
run_options["T"] = args.T
run_options["h"] = args.period
if args.backend == "vnative":
    import nonos.bridge

    vnative_options = nonos.bridge.SpaceFilterOptions()
    if args.neighborhood_radius is not None:
        vnative_options.neighborhood_radius = args.neighborhood_radius
    if args.min_radius is not None:
        vnative_options.min_radius = args.min_radius
    run_options["vnative_options"] = vnative_options

# Run the simulation from the inputs previously defined and add
# results to the hdf5 file. The visualisation of the output may be done
# with the vview command.
with MechanicsHdf5Runner(
    mode="r+", io_filename=args.file[0], config=RunnerConfig(args.backend)
) as io:
    # By default earth gravity is applied and the units are those
    # of the International System of Units.
    # Because of fixed collision margins used in the collision detection,
    # sizes of small objects may need to be expressed in cm or mm.
    io.run(run_options)
