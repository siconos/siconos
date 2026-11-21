#!/usr/bin/env @Python_EXECUTABLE@
"""
Description: SVG pictures of the output steps of a 2D Siconos mechanics-IO
HDF5 file: the static objects, the fixed segments hatched on their material
side, the bodies coloured by their speed, the contacts (reaction impulse,
Coulomb cone).  One file per picture, <file>_svg_<time>.svg beside the HDF5
file.  --time draws the output step at or before the time; --frame an
output step index (negative from the end); the last step by default.

Simulations run with run_options["svg_output"] = "0.5,1.0" (or --svg-output
on the siconos command line) write the same pictures during the run.

PNG: rsvg-convert -o picture.png picture.svg (guix shell librsvg).
"""

import argparse

from siconos.io.svg import convert


def material_side(text):
    """ID=+1 or ID=-1"""
    key, _, value = text.partition("=")
    if value not in ("+1", "1", "-1"):
        raise argparse.ArgumentTypeError(f"{text}: ID=+1 or ID=-1")
    return int(key), int(value)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter
    )
    parser.add_argument("io_filename", help="the siconos HDF5 file")
    parser.add_argument(
        "--time", type=float, action="append", default=[], help="a time, repeatable"
    )
    parser.add_argument(
        "--frame",
        type=int,
        action="append",
        default=[],
        help="an output step index, repeatable",
    )
    parser.add_argument("-o", "--output", help="the SVG file (one picture only)")
    parser.add_argument("--dir", help="the directory of the pictures")
    parser.add_argument(
        "--mode",
        choices=("auto", "circle", "use", "bins"),
        default="auto",
        help="one element per body (circle), one <use> of its shape (use), the "
        "bodies merged by colour (bins); auto by the number of bodies",
    )
    for name, what in (
        ("arrows", "the reaction arrows"),
        ("cones", "the Coulomb cones"),
    ):
        group = parser.add_mutually_exclusive_group()
        group.add_argument(
            f"--{name}",
            dest=name,
            action="store_true",
            default=None,
            help=f"draw {what} (default: by the number of contacts)",
        )
        group.add_argument(f"--no-{name}", dest=name, action="store_false")
    parser.add_argument(
        "--material-side",
        type=material_side,
        action="append",
        default=[],
        metavar="ID=+1|-1",
        help="the material side of the segments of the static object ID: on the "
        "side of perp(p2 - p1) (+1) or the other (-1), instead of the inference",
    )
    parser.add_argument("--labels", action="store_true", help="the body ids")
    args = parser.parse_args()

    names = convert(
        args.io_filename,
        times=args.time,
        frames=args.frame,
        directory=args.dir,
        output=args.output,
        mode=args.mode,
        arrows=args.arrows,
        cones=args.cones,
        material_sides=dict(args.material_side),
        labels=args.labels,
    )
    for name in names:
        print(name)
