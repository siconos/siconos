#!/usr/bin/env @Python_EXECUTABLE@
"""
Description: Show information about a Siconos mechanics-IO HDF5 file.
"""

# Lighter imports before command line parsing
import argparse

parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument(
    "file", metavar="input", type=str, nargs="+", help="input file(s) (HDF5)"
)
parser.add_argument(
    "-O", "--list-objects", action="store_true", help="List object names in the file"
)
parser.add_argument(
    "-C",
    "--list-contactors",
    action="store_true",
    help="List contactor names in the file",
)
parser.add_argument("-V", "--version", action="version", version="@SICONOS_VERSION@")

if __name__ == "__main__":
    args = parser.parse_args()

# Heavier imports after command line parsing
import numpy as np
from siconos.io.mechanics_hdf5 import MechanicsHdf5


def summarize(io):
    io.static_data()
    dpos_data = io.dynamic_data()
    cf_data = io.contact_forces_data()
    solv_data = io.solver_data()
    t0 = dpos_data[:, 0].min()
    t1 = dpos_data[:, 0].max()
    times, counts = np.unique(dpos_data[:, 0], return_counts=True)
    print(f"Time simulated: {t0} to {t1} = {len(times)} steps")

    cf_times, cf_counts = np.unique(cf_data[:, 0], return_counts=True)
    min_cf = 0
    if len(cf_counts) != 0:
        min_cf = cf_counts.min()

    # Are there times where there are no contact forces?
    if len(np.setdiff1d(times, cf_times, assume_unique=True)) > 0:
        min_cf = 0

    print()
    print("            {:>10} {:>10} {:>10}".format("Min", "Avg", "Max"))
    print("            {:->10} {:->10} {:->10}".format("", "", ""))
    print(
        f"Objects:    {counts.min(): >10} {int(counts.mean()): >10} {counts.max(): >10}"
    )
    if len(cf_counts) != 0:
        print(
            f"Contacts:   {min_cf: >10} {int(cf_counts.mean()): >10} {cf_counts.max(): >10}"
        )
    else:
        print(f"Contacts:   {min_cf: >10} {0: >10} {0: >10}")

    print(
        f"Iterations: {int(solv_data[:, 1].min()): >10} {int(solv_data[:, 1].mean()): >10} {int(solv_data[:, 1].max()): >10}"
    )
    print(
        f"Precision:  {solv_data[:, 2].min(): >10.3g} {solv_data[:, 2].mean(): >10.3g} {solv_data[:, 2].max(): >10.3g}"
    )
    print(
        f"Loc. Prec.: {solv_data[:, 3].min(): >10.3g} {solv_data[:, 3].mean(): >10.3g} {solv_data[:, 3].max(): >10.3g}"
    )


def list_objects(io):
    print()
    print("Objects:")
    print()
    print("{:>5} {:>15} {:>6} {:>6}".format("Id", "Name", "Mass", "ToB"))
    print("{0:->5} {0:->15} {0:->6} {0:->6}".format(""))
    for name, obj in io.instances().items():
        print(
            "{:>5} {:>15} {:>6.4g} {:>6.4g}".format(
                obj.attrs["id"], name, obj.attrs["mass"], obj.attrs["time_of_birth"]
            )
        )


def list_contactors(io):
    print()
    print("Contactors:")
    print()
    print("{:>5} {:>15} {:>9} {:>9}".format("Id", "Name", "Type", "Primitive"))
    print("{0:->5} {0:->15} {0:->9} {0:->9}".format(""))
    for name, obj in io.shapes().items():
        print(
            "{:>5} {:>15} {:>9} {:>9}".format(
                obj.attrs["id"],
                name,
                obj.attrs["type"],
                obj.attrs.get("primitive", ""),
            )
        )


def compute_violation(io):
    cf_data = io.contact_forces_data()
    gap = cf_data[:, 14]
    negative_gap = np.where(gap < 0, gap, 0.0)

    print("            {:>10} {:>10} {:>10}".format("Min", "Avg", "std"))

    print(
        f"Violation:  {negative_gap.min(): >10.2e} {negative_gap.mean(): >10.2e} {negative_gap.std(): >10.2e}"
    )


if __name__ == "__main__":
    try:
        with MechanicsHdf5(mode="r", io_filename=args.file[0]) as io:
            if io.dynamic_data() is None or len(io.dynamic_data()) == 0:
                print("Empty simulation found.")
            else:
                print()
                print(f'Filename: "{args.file[0]}"')
                summarize(io)
                compute_violation(io)
                if args.list_objects:
                    list_objects(io)
                if args.list_contactors:
                    list_contactors(io)
    except OSError as e:
        print(f'Error reading "{args.file[0]}"')
        print(e)
