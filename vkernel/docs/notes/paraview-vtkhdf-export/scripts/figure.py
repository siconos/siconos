# The figures of paraview-vtkhdf-export.md: a screenshot of a ParaView state
# written with the VTKHDF files, at a time, with or without the contacts.
#
#   pvpython --force-offscreen-rendering figure.py <name>_paraview.py out.png
#       [--time T] [--contacts] [--bounds xmin xmax ymin ymax]
#       [--size width height]

import argparse

parser = argparse.ArgumentParser()
parser.add_argument("state")
parser.add_argument("output")
parser.add_argument("--time", type=float)
parser.add_argument("--contacts", action="store_true")
parser.add_argument("--bounds", type=float, nargs=4)
parser.add_argument("--size", type=int, nargs=2, default=(1600, 700))
options = parser.parse_args()

# the state defines bodies, view, display, and contacts when written with them
exec(open(options.state).read())

view.ViewSize = list(options.size)
if options.time is not None:
    GetAnimationScene().AnimationTime = options.time
if options.contacts:
    Show(contacts, view)
    # the range grows with the steps shown: at one step, its own range
    Render(view)
    contacts_display.RescaleTransferFunctionToDataRange(False, True)
    contacts_display.SetScalarBarVisibility(view, True)
    display.Opacity = 0.3
if options.bounds:
    view.ResetCamera(*options.bounds, -1.0, 1.0)
Render(view)
SaveScreenshot(options.output, view)
