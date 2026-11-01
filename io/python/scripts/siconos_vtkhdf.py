#!/usr/bin/env @Python_EXECUTABLE@
"""
Description: Write the bodies of a Siconos mechanics-IO HDF5 file in a
transient VTKHDF file, readable by ParaView, with the shapes of the bodies
(<output>_shapes.vtkhdf), the static objects (<output>_static.vtkhdf), the
contacts (<output>_contacts.vtkhdf), the finite element meshes
(<output>_fem.vtkhdf) and a ParaView state drawing them
(paraview --state=<output>_paraview.py).

Simulations run with run_options["output_format"] = "vtkhdf" (or "both")
write this file directly.
"""

import argparse
import os

from siconos.io.vtkhdf import convert

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("io_filename", help="the siconos HDF5 file")
    parser.add_argument(
        "-o",
        "--output",
        help="the VTKHDF file, <io_filename without extension>.vtkhdf by default",
    )
    parser.add_argument(
        "--compression", action="store_true", help="gzip the VTKHDF datasets"
    )
    args = parser.parse_args()

    output = args.output or f"{os.path.splitext(args.io_filename)[0]}.vtkhdf"
    steps = convert(args.io_filename, output, compression=args.compression)
    print(f"{output}: {steps} steps")
