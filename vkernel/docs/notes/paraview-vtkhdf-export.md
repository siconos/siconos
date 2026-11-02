# Viewing vkernel simulations in ParaView: the VTKHDF export

How to write the results of a mechanics simulation (vnative backend, and
the other backends through the same runner) in files that ParaView reads
natively, and how to open them. The export is made for large granular
simulations, 100k to 1M bodies, where vview and the old VTK export
(`siconos_vexport`) do not scale: a body is one point, drawn by the GPU
with its shape, and a time step is a contiguous block of an HDF5 file.

Code: `io/python/vtkhdf.py` (the writers, the shapes, the converter),
`io/python/mechanics_run.py` (`output_format`, `output_bodies`,
`output_contact_forces`, `output_static_vtkhdf`, `output_fem_vtkhdf`),
`io/python/scripts/siconos_vtkhdf.py` (the command). Introduced in
63018115d.

![the dam break of bench-dambreak.py 40 at t = 0.9: disks and rocks drawn
with their shape, colored by velocity, the terrain in gray](paraview-vtkhdf-export/dambreak_flow.png)

## 1. Export

### During the simulation

One run option, before `io.run(run_options)`:

```python
run_options["output_format"] = "vtkhdf"  # or "both", default "siconos"
```

| `output_format`       | what is written                                                                                         |
|-----------------------|---------------------------------------------------------------------------------------------------------|
| `"siconos"` (default) | the HDF5 file of the simulation, unchanged                                                              |
| `"vtkhdf"`            | the ParaView files **instead of** the tables `dynamic`, `velocities`, `radii` and `cf` of the HDF5 file |
| `"both"`              | the HDF5 file as by default, and the ParaView files                                                     |

With `"vtkhdf"` the HDF5 file still holds the scene (`input`, `ref`,
nonsmooth laws), the static objects, the FEM tables, the contact info
and work, the logs: vview cannot show the run, but a restart works (§4).

Without editing the script:

```
python3 bench-dambreak.py 40 --output-format vtkhdf
siconos --output-format vtkhdf bench-disks.py vnative 50
```

`bench-dambreak.py` has the option; the `siconos` command
(`--t0 --T --h --output-format`, and the vnative time stepping options)
sets them for any python driver it runs, over what the driver chose.

### After the simulation

Any siconos HDF5 file, old ones included:

```
siconos_vtkhdf my_run.hdf5                 # my_run.vtkhdf and the others
siconos_vtkhdf my_run.hdf5 -o out.vtkhdf   # out.vtkhdf, out_shapes.vtkhdf...
siconos_vtkhdf my_run.hdf5 --compression   # gzip the datasets
```

In the guix container (`guix shell --container`) the installed scripts
do not start, there is no `/usr/bin/env`: `python3 <install>/bin/siconos_vtkhdf ...`.

## 2. The files

Next to `my_run.hdf5` (or the `-o` name):

| file                     | content                                       | written                                              |
|--------------------------|-----------------------------------------------|------------------------------------------------------|
| `my_run.vtkhdf`          | the dynamic bodies, one point per body        | always                                               |
| `my_run_shapes.vtkhdf`   | the glyph table: one block per distinct shape | always                                               |
| `my_run_static.vtkhdf`   | the static objects alive at each step         | if there are static objects                          |
| `my_run_contacts.vtkhdf` | the contacts, one point per contact           | if the contacts are output (`output_contact_forces`) |
| `my_run_fem.vtkhdf`      | the finite element meshes, deformed           | if there are FEM bodies                              |
| `my_run_paraview.py`     | a ParaView state drawing all of them          | always                                               |

All the files have the same time steps: a step without contact is an
empty step of the contacts file, so ParaView never shows the contacts of
another time.

The point and cell arrays:

| file     | arrays                                                                                                                                                           |
|----------|------------------------------------------------------------------------------------------------------------------------------------------------------------------|
| bodies   | `id`, `orientation` (quaternion w, x, y, z), `velocity`, `angular_velocity`, `radius` (vnative: disk radius, polygon circumradius), `shape_index`, `glyph_scale` |
| static   | `id`, `orientation`, `shape_index`, `glyph_scale`                                                                                                                |
| contacts | `normal`, `reaction` (impulse, global frame), `reaction_normal`, `reaction_tangent` (norm), `gap`, `relative_velocity`, `mu`, `interaction`, `ds1`, `ds2`        |
| FEM      | points: `displacement`; cells: `strain`, `stress` (xx, yy, xy), `von_mises` (plane stress), `body` (ds id)                                                       |

A contact is drawn at the middle of its two contact points. In 2D the
z coordinates are 0 and the orientation is a rotation about z.

## 3. View

```
paraview --state=my_run_paraview.py
```

or *File > Load State* in ParaView, or `pvpython my_run_paraview.py`.
The state:

- draws the bodies with the *3D Glyphs* representation: each body takes
  the block `shape_index` of `my_run_shapes.vtkhdf`, rotated by
  `orientation` and scaled by `glyph_scale` (the radius for the disks and
  spheres, which share a unit block; 1 for the others), colored by the
  velocity;
- draws the static objects the same way, in gray;
- draws the FEM meshes with their edges, colored by the von Mises stress;
- loads the contacts **hidden**: the eye of `contacts` in the Pipeline
  Browser shows them, as small spheres colored by the reaction;
- frames the camera on the static objects and the initial positions of
  the bodies.

The color ranges grow with the steps shown. The state holds absolute
paths: after moving the files, convert again or edit the paths at the top
of the state.

![the contacts at t = 0.45, the bodies made transparent, colored by the
reaction impulse](paraview-vtkhdf-export/dambreak_contacts.png)

Other ways to look at the same files, in ParaView:

- disks or spheres only: the *Point Gaussian* representation of the
  bodies, *Scale By Array* `radius`, *Use Scale Function* off, is lighter
  than the glyphs;
- the force chains: color the contacts by `reaction_normal`, or a
  *Threshold* on it;
- the bodies of a contact: `ds1`, `ds2`; a body: *Find Data* on `id`.

## 4. How it works, and what it costs

**One point per body.** The bodies file is a transient VTKHDF PolyData:
every step appends its rows (positions, then each array) to the same
datasets, and `Steps` records where each step starts. The vertex
topology, one vertex per point, is written once as `0..max(N)` and every
step reads a prefix of it, so the number of bodies can change between
steps (births, deaths) at no cost. The writer needs numpy and h5py
only, not VTK.

**The glyph table.** `ShapeTable` reads `input` and `ref`: an
`add_object` body takes the block of its shape, or of its contactors
merged when there are several or one is offset; the grain k of an
`add_objects` aggregate has the id of the aggregate + k (as vview), and
the block of its shape. Shapes: 2D `Disk`, `Circle`, `Box2d`, convex
polygons, `Segment`, `Line` (drawn 1000 long); 3D `Sphere`, `Box`,
`Cylinder`, `Cone`, `Capsule` (as vview draws them), convex hulls. Any
other shape is drawn as a small cross, with a warning. The convex
vertices are taken as given: vnative requires them centred on their
centroid, so they are in the body frame.

**Static objects.** One point per static object at each step, from the
`static` table, kept only from its `time_of_birth` to its
`time_of_death` (the table keeps the dead ones).

**FEM.** The objects with a `material` are FEM bodies (stored with
type `static`). Their mesh is the gmsh 2.2 text of their shape; its
points are put in the "spatial order" of `mechanics_run` (vertices
sorted by (y, x)), the order of the rows of `fem_displacements_<id>`.
The triangles (T3) are written once; each step writes the deformed
points and the arrays. The FEM tables have no row at the initial
output: the first step of the file is the undeformed mesh.

**Step index.** Since 63018115d the HDF5 file has `data/step_index`:
for each table that gets rows at each output step, one row
`(time, first row, number of rows)` per step, empty steps included.
vview, `radii(step)`, `p0s(step)`, the cold restart and the converter
use it, and fall back on the time column for the files without it. A
file whose tables were written without an index does not get one when
a run is appended to it.

**Restart.** A run continued from its HDF5 file (cold restart) appends
to the VTKHDF files; in `"vtkhdf"` mode the last positions and
velocities come from `my_run.vtkhdf`. The grains of `add_objects`
aggregates restart from the last output too.

**Costs** (measured, numpy writer, one core):

| bodies | write / step | disk / step | ParaView load of a step |
|--------|--------------|-------------|-------------------------|
| 100k   | 21 ms        | 6.2 MB      | 4 ms                    |
| 1M     | 153 ms       | 62 MB       | 21 ms                   |

The shape table of 1.1M bodies (102 shapes) builds in 0.06 s and adds
23 ms per step. The positions are float32, as the siconos tables.
`bench-dambreak.py 40` (1136 bodies, 2985 steps, 101 outputs) runs in
2 min 16 s and writes 9.6 MB of bodies and 9.7 MB of contacts.

## 5. Limits

- FEM: 2D T3 meshes only; no T6, quadrangles, 3D elements.
- The strain and stress of an element are matched to the triangles in
  the order of the gmsh file, as `computeStrainTensor` produces them; the
  counts are checked, the order is assumed.
- Mesh shapes of rigid bodies (`vtp`, `stl`, `brep` files) and heightmaps
  are drawn as crosses.
- The FEM tables are written in every mode.
- `siconos_vexport` (`vexport.py`) is still there, outdated and slow.

## 6. Reproduce the figures

```
python3 bench-dambreak.py 40 --t-settle 0.3 --T 1.0 --output-format vtkhdf
pvpython --force-offscreen-rendering figure.py bench-dambreak-40_paraview.py \
    dambreak_flow.png --time 0.9 --bounds 0 700 60 560 --size 1600 900
pvpython --force-offscreen-rendering figure.py bench-dambreak-40_paraview.py \
    dambreak_contacts.png --time 0.45 --contacts --bounds 40 130 320 400 \
    --size 1200 900
```

`paraview-vtkhdf-export/scripts/figure.py`: a screenshot of a state at a
time, with the contacts or not. ParaView 6.0 from guix
(`guix shell paraview`).
