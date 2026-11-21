#!/usr/bin/env @Python_EXECUTABLE@
"""
Description: SVG pictures of the output steps of a 2D Siconos mechanics-IO
HDF5 file: the static objects, the fixed segments hatched on their material
side, the bodies coloured by their speed, the contacts (reaction impulse,
Coulomb cone, force chains), the FEM meshes deformed.  One file per
picture, <file>_svg_<time>.svg beside the HDF5 file.  --time draws the
output step at or before the time; --frame an output step index (negative
from the end); the last step by default.  --coulomb adds the (r_n, r_t)
plot of each picture, <file>_coulomb_<time>.svg.

--animate writes <file>_svg_animation.svg, the output steps animated
(SMIL: a browser plays it, a viewer without SMIL shows the first step);
--strobe <file>_svg_strobe.svg, the steps superimposed.  With either, the
pictures of --time and --frame are drawn in the view of the animation:
--frame N is the frame N of the animation (of the steps kept by --stride).

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
        ("chains", "the force chains"),
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
    parser.add_argument(
        "--scale",
        choices=("proportional", "fixed", "log"),
        default="proportional",
        help="the length of the arrows: proportional to the impulse, fixed (the "
        "direction only) or log (3 decades)",
    )
    parser.add_argument(
        "--colour",
        choices=("speed", "omega", "none"),
        default="speed",
        help="the colouring of the bodies: |v|, |omega| or none",
    )
    parser.add_argument("--labels", action="store_true", help="the body ids")
    parser.add_argument(
        "--view",
        type=float,
        nargs=4,
        metavar=("X0", "X1", "Y0", "Y1"),
        help="crop the pictures to this window (said in their title)",
    )
    parser.add_argument(
        "--coulomb", action="store_true", help="the Coulomb plot of each picture"
    )
    parser.add_argument(
        "--animate", action="store_true", help="the animation of the output steps"
    )
    parser.add_argument(
        "--strobe", action="store_true", help="the output steps superimposed"
    )
    parser.add_argument(
        "--stride",
        type=int,
        default=1,
        help="every stride-th output step (and the last)",
    )
    parser.add_argument(
        "--uniform",
        action="store_true",
        help="one frame per output step, instead of the physical times",
    )
    parser.add_argument(
        "--discrete", action="store_true", help="no interpolation between the steps"
    )
    parser.add_argument(
        "--speed", type=float, default=1.0, help="the playback speed, 1: real time"
    )
    parser.add_argument(
        "--controls",
        action="store_true",
        help="buttons under the animation: pause, step, seek (a script)",
    )
    args = parser.parse_args()

    names = convert(
        args.io_filename,
        times=args.time,
        frames=args.frame,
        directory=args.dir,
        output=args.output,
        mode=args.mode,
        coulomb=args.coulomb,
        animate=args.animate,
        strobe=args.strobe,
        stride=args.stride,
        uniform=args.uniform,
        discrete=args.discrete,
        speed=args.speed,
        controls=args.controls,
        arrows=args.arrows,
        chains=args.chains,
        cones=args.cones,
        scale=args.scale,
        colour=args.colour,
        material_sides=dict(args.material_side),
        labels=args.labels,
        view=args.view,
    )
    for name in names:
        print(name)
