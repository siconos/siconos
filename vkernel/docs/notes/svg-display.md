# SVG pictures of a 2D simulation

A still picture of a 2D mechanics simulation at chosen times: the static
objects, the bodies with their shape coloured by their speed, the contacts
(the reaction impulse, the Coulomb cone), and the fixed segments hatched on
their material side, so that a floor reads as a floor. One self-contained
SVG file per picture, for figures and for checking a scene by eye; for
large runs and animations, use the ParaView export
([paraview-vtkhdf-export.md](paraview-vtkhdf-export.md)).

Code: `io/python/svg.py` (the picture, the side inference, the converter,
the pictures of a run), `io/python/mechanics_run.py` (`svg_output`,
`output_svg`), `io/python/scripts/siconos_svg.py` (the command). Test:
`vkernel/src/siconos/vkernel/tests/test_svg.py` (ctest
`vkernel_python_svg`).

![two disks rolling down a slope against a wall, t = 1: the slope hatched
below, the wall on its right; at rest against the wall, the second disk
with its two contacts](svg-display/slope_wall.svg)

## 1. Making the pictures

### During the simulation

```python
run_options["svg_output"] = "0.5,1.0,120"
```

or, for any python driver, `siconos --svg-output 0.5,1.0,120 driver.py`.
A value with a decimal point (or an exponent) is a time, an int a
simulation step number (its time is k h). The picture is made at the step
k with t_k <= request < t_k + h, whether it is an output step or not; a
request the run does not reach is reported on stderr. The files are
`<hdf5 file>_svg_<time>.svg` beside the HDF5 file, the time with 6
decimals and zero padded (`run_svg_000001.250000.svg`) so that they sort.

The run draws from the arrays it writes in the HDF5 tables, rounded to
float32 as the tables are: it works with `output_format` `"vtkhdf"` (no
`dynamic` table), and its picture of an output step is the same bytes as
the converter's picture of that step.

### After the simulation

```
siconos_svg run.hdf5                          # the last output step
siconos_svg run.hdf5 --time 0.5 --time 1.0    # the output steps at or before
siconos_svg run.hdf5 --frame 0 --frame -1     # output step indices
siconos_svg run.hdf5 --time 0.27 -o fig.svg
siconos_svg run.hdf5 --dir pictures
```

A file only has its output steps: a time between two of them draws the
one before, never an interpolation, and the title gives both times
(`t = 0.25 s (requested 0.27)`); a time before the first output step is
an error. The file needs a step index (`data/step_index`, the files of
2026 on).

In the guix container the installed scripts do not start (no
`/usr/bin/env`): `python3 <install>/bin/siconos_svg run.hdf5`.

### PNG

A conversion, not an output: `rsvg-convert -o fig.png fig.svg` (in the
guix container: `guix shell librsvg -- rsvg-convert ...`), or Inkscape,
or cairosvg.

## 2. What is drawn

In this order, each layer in a fixed order so that a file is the same
bytes from one run to the next:

1. the static objects alive at the time (`static_geometry` of
   `vtkhdf.ShapeTable`, as in ParaView): closed polygons filled grey,
   segments as lines, hatched on their material side (section 3);
2. the bodies by id, with their shape (the blocks of `ShapeTable`: disks,
   convex polygons, Box2d, composite bodies), rotated by their angle,
   coloured by |v| on a fixed five-colour ramp; the legend gives the range
   of the drawn step;
3. the contacts by (ds1, ds2, interaction): the Coulomb cone (half-angle
   atan(mu) about the normal, at the point of ds1) and the reaction
   impulse on ds1 (column `reaction` of `cf`, H^T lambda), scaled so that
   the largest of the step is 1.2 typical body sizes; the title gives
   that largest impulse in N s. It is an impulse under Moreau-Jean, not a
   force.

The title gives the time drawn (and the one requested), the integrator
(`osi` of the run options, written in
`data/siconos_mechanics_run_options` since the run options are output
again: files of older runs have none), the numbers of bodies and contacts
and the friction coefficients.

| bodies        | mode (`--mode auto`)                                         |
|---------------|--------------------------------------------------------------|
| up to 2000    | `circle`: one element per body, a rotation tick on the disks |
| up to 20000   | `use`: each shape once in `<defs>`, one `<use>` per body     |
| above         | `bins`: the bodies merged in one path per colour, no outline |

The arrows are drawn up to 2000 contacts, the cones up to 20 (they are
illegible beyond); `--arrows`, `--no-arrows`, `--cones`, `--no-cones`
force them. `--labels` writes the body ids.

Measured on a lattice of disks (bench_scale.py, siconos.asor-work/svg;
converter time on the CPU, rsvg-convert at 1600 px):

| bodies | circle, with arrows      | use                       | bins                      |
|--------|--------------------------|---------------------------|---------------------------|
| 1000   | 0.33 MB, 0.05 s          | 0.055 MB, 0.009 s         | 0.063 MB, 0.010 s         |
| 10000  | 3.3 MB, 0.28 s, rsvg 1.1 s | 0.55 MB, 0.04 s, rsvg 0.65 s | 0.60 MB, 0.05 s, rsvg 0.54 s |
| 100000 | 33 MB, 3.4 s, rsvg 6.1 s | 5.5 MB, 0.37 s, rsvg 3.2 s | 6.0 MB, 0.51 s, rsvg 2.0 s |

The merged paths are not smaller than the `<use>` (each arc repeats the
radius), they render faster.

## 3. The material side of a segment

A segment has no side in the model (`collision/shape/segment.hpp`). The
side is inferred:

- from the contacts with the segment. A contact with a fixed obstacle has
  ds1 == ds2 (the body: the obstacle is not named in `cf`) and its
  `point_a` on the obstacle; it is matched to the segment by distance
  (1e-5 of the scene diagonal). Its normal points from the obstacle into
  the body: the material is on the other side. A contact at an end point,
  its normal nearly along the segment, does not count. The contacts are
  those of the output steps up to the drawn one, so that a hatch does not
  change with a contact that comes and goes;
- without contact, from the bodies near the segment (centre within two
  typical sizes, projection inside the segment), at the drawn step then
  at the first output step;
- with bodies on both sides (a wall between two compartments), or none,
  the segment is not hatched and a warning names it (static id and end
  points).

When the inference is wrong (a ledge hit once from below is hatched on
top), give the side: `--material-side ID=+1` hatches the segments of the
static object ID (a negative id) on the side of perp(p2 - p1), p1 and p2
the end points of the segment, `-1` on the other. A side declared in the
scene would be a change of the HDF5 format, left to the hdf5_v2 work.
