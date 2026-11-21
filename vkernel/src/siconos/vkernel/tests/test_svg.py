# Siconos is a program dedicated to modeling, simulation and control
# of non smooth dynamical systems.
#
# Copyright 2026 INRIA.
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
# http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.
#

"""The SVG pictures of siconos.io.svg on vnative runs.

slope: two disks roll down a slope against a wall.
- the pictures of the run (svg_output, svg_options) and of the file
  (convert) are the same bytes, the animation and the Coulomb plot too, with
  output_format "vtkhdf" (no dynamic tables) as well;
- the slope is hatched below, the wall on its right, from the contacts;
- a disk is drawn where the dynamic table puts it, its contact points on
  its outline;
- a time between two output steps draws the step before it;
- the picture at t = 1 is the reference svg_slope_wall.ref.svg
  (SICONOS_SVG_UPDATE_REF=1 rewrites it).
wall: a disk thrown at a wall from each side: the wall is not hatched.
bounce: a disk thrown on a floor, animated: the arrow tip is on an edge of
the cone exactly when the contact slides; the picture of a frame and the
frame of the animation agree; the hatch is the same in all the frames.
restart: the slope run to 0.5 then continued to 1: the converter draws the
picture of the continued run.
The side inference on synthetic contacts: both sides, end points, the
override.

Each simulation runs in its own process (a second vnative simulation in a
process crashes). Run with pytest, nonos and siconos.io in the PYTHONPATH.
"""

import os
import re
import subprocess
import sys
import xml.etree.ElementTree as ET

import h5py
import numpy as np
import pytest

from siconos.io import svg

HERE = os.path.dirname(os.path.abspath(__file__))
REFERENCE = os.path.join(HERE, "svg_slope_wall.ref.svg")
RADIUS = 0.5
TIMES = "0.5,1.0,0.27"
MU_BOUNCE = 0.1


def build(io, scene):
    from siconos.mechanics.collision.tools import Contactor

    io.add_primitive_shape("Disk", "Disk", [RADIUS])
    if scene == "slope":
        io.add_primitive_shape("Slope", "Segment", (-3, 1.5, 3, 0))
        io.add_primitive_shape("Wall", "Segment", (3, 0, 3, 2))
        io.add_object("slope", [Contactor("Slope")], translation=[0, 0])
        io.add_object("wall", [Contactor("Wall")], translation=[0, 0])
        # on the slope y = 0.75 - x / 4, a little above it
        disks = [((x, 0.75 - 0.25 * x + 0.55), (0, 0)) for x in (0.0, 1.8)]
        io.add_Newton_impact_friction_nsl("contact", mu=0.3, e=0)
    elif scene == "wall":
        io.add_primitive_shape("Floor", "Segment", (-3, 0, 3, 0))
        io.add_primitive_shape("Wall", "Segment", (0, 0, 0, 2))
        io.add_object("floor", [Contactor("Floor")], translation=[0, 0])
        io.add_object("wall", [Contactor("Wall")], translation=[0, 0])
        disks = [((-1.5, RADIUS), (3, 0)), ((1.5, RADIUS), (-3, 0))]
        io.add_Newton_impact_friction_nsl("contact", mu=0.3, e=0)
    else:
        io.add_primitive_shape("Floor", "Segment", (-1, 0, 16, 0))
        io.add_object("floor", [Contactor("Floor")], translation=[0, 0])
        disks = [((0, 1.5), (6, 0))]
        io.add_Newton_impact_friction_nsl("contact", mu=MU_BOUNCE, e=0.5)
    for k, (x, v) in enumerate(disks):
        io.add_object(
            f"disk{k}",
            [Contactor("Disk")],
            translation=list(x),
            orientation=[0],
            velocity=[v[0], v[1], 0],
            mass=1,
            inertia=0.5 * RADIUS**2,
        )


def run(filename, output_format, scene, T=1.0, new=True):
    """A scene, run with the pictures of TIMES, its animation and Coulomb
    plots; new: build the scene, else continue the run of the file."""
    import siconos.numerics as sn
    from siconos.io.mechanics_run import (
        MechanicsHdf5Runner,
        MechanicsHdf5Runner_run_options,
        RunnerConfig,
    )
    import nonos  # noqa: F401 registers the vnative backend

    config = RunnerConfig("vnative")
    if new:
        with MechanicsHdf5Runner(config=config, io_filename=filename) as io:
            build(io, scene)

    options = sn.SolverOptions(sn.solver_ids.SICONOS_FRICTION_2D_NSGS)
    options.iparam[sn.params.SICONOS_IPARAM_MAX_ITER] = 1000
    options.dparam[sn.params.SICONOS_DPARAM_TOL] = 1e-10
    run_options = MechanicsHdf5Runner_run_options()
    run_options["t0"] = 0
    run_options["T"] = T
    run_options["h"] = 0.005
    run_options["solver_options"] = options
    run_options["output_frequency"] = 1 if scene == "bounce" else 10
    run_options["output_format"] = output_format
    run_options["svg_output"] = TIMES
    run_options["svg_options"] = {"animate": True, "coulomb": True}
    run_options["verbose"] = False
    with MechanicsHdf5Runner(mode="r+", config=config, io_filename=filename) as io:
        io.run(run_options)


def simulate(directory, output_format="siconos", scene="slope", T=1.0, new=True):
    os.makedirs(directory, exist_ok=True)
    filename = os.path.join(directory, f"{scene}.hdf5")
    args = [filename, output_format, scene, str(T), "new" if new else "continue"]
    subprocess.run([sys.executable, __file__, *args], check=True, cwd=directory)
    return filename


@pytest.fixture(scope="module")
def siconos_run(tmp_path_factory):
    return simulate(str(tmp_path_factory.mktemp("siconos")))


@pytest.fixture(scope="module")
def bounce_run(tmp_path_factory):
    # it bounces and slides, then rolls from about t = 2.3
    return simulate(str(tmp_path_factory.mktemp("bounce")), scene="bounce", T=2.6)


def read(name):
    with open(name) as f:
        return f.read()


def picture(filename, time):
    with open(svg.svg_filename(filename, time)) as f:
        return f.read()


def test_run_and_file_same_bytes(siconos_run, tmp_path):
    for time in (0.5, 1.0):
        (name,) = svg.convert(siconos_run, times=[time], directory=str(tmp_path))
        with open(name) as f:
            assert f.read() == picture(siconos_run, time)
    # 0.27 is not an output step: drawn by the run at its own step
    assert "t = 0.27 s" in picture(siconos_run, 0.27)


def test_deterministic(siconos_run, tmp_path):
    first = svg.convert(siconos_run, output=str(tmp_path / "a.svg"))[0]
    second = svg.convert(siconos_run, output=str(tmp_path / "b.svg"))[0]
    with open(first) as a, open(second) as b:
        assert a.read() == b.read()


def _hatch(text, static_id):
    """The first stroke (dx, dy in px, y down) of the hatch of a static
    object, and its line."""
    path = re.search(
        rf'<path id="h{-static_id}_0" d="M[^l]*l([-0-9.]+),([-0-9.]+)', text
    )
    line = re.search(
        rf'<line id="l{-static_id}_0" x1="([-0-9.]+)" y1="([-0-9.]+)" '
        rf'x2="([-0-9.]+)" y2="([-0-9.]+)"',
        text,
    )
    assert path is not None and line is not None
    return np.array(path.groups(), float), np.array(line.groups(), float)


def test_material_sides(siconos_run):
    text = picture(siconos_run, 1.0)
    with h5py.File(siconos_run, "r") as f:
        ids = {name: int(o.attrs["id"]) for name, o in f["data/input"].items()}
    # the slope: hatched below (y down in the picture)
    stroke, _ = _hatch(text, ids["slope"])
    assert stroke[1] > 0
    # the wall at x = 3, the disks on its left: hatched on the right
    stroke, _ = _hatch(text, ids["wall"])
    assert stroke[0] > 0


def test_body_pose(siconos_run):
    """The disks where the table puts them: the slope end points give the
    affine map of the picture."""
    text = picture(siconos_run, 1.0)
    with h5py.File(siconos_run, "r") as f:
        data = f["data"]
        ids = {name: int(o.attrs["id"]) for name, o in data["input"].items()}
        steps = svg.FileSteps(data)
        bodies = steps.step(steps.resolve(time=1.0)).bodies
    _, (x1, y1, x2, y2) = _hatch(text, ids["slope"])
    scale = (x2 - x1) / 6.0
    # y down in the picture
    assert (y2 - y1) / 1.5 == pytest.approx(scale, rel=1e-3)
    for row in bodies:
        circle = re.search(
            rf'<g id="b{int(row[0])}" fill="[^"]*"><circle cx="([-0-9.]+)" '
            rf'cy="([-0-9.]+)" r="([-0-9.]+)"',
            text,
        )
        cx, cy, r = map(float, circle.groups())
        assert cx == pytest.approx(x1 + (row[1] + 3) * scale, abs=0.02)
        assert cy == pytest.approx(y1 - (row[2] - 1.5) * scale, abs=0.02)
        assert r == pytest.approx(RADIUS * scale, abs=0.01)


def test_time_between_output_steps(siconos_run, tmp_path):
    (name,) = svg.convert(siconos_run, times=[0.27], directory=str(tmp_path))
    assert name.endswith("_svg_000000.250000.svg")
    with open(name) as f:
        assert "t = 0.25 s (requested 0.27)" in f.read()
    (name,) = svg.convert(siconos_run, frames=[-1], directory=str(tmp_path))
    assert name.endswith("_svg_000001.000000.svg")
    with pytest.raises(ValueError):
        svg.convert(siconos_run, times=[-1.0], directory=str(tmp_path))


def test_reference(siconos_run):
    text = picture(siconos_run, 1.0)
    if os.environ.get("SICONOS_SVG_UPDATE_REF"):
        with open(REFERENCE, "w") as f:
            f.write(text)
    with open(REFERENCE) as f:
        assert text == f.read()


def test_vtkhdf_output_format(siconos_run, tmp_path):
    """No dynamic table: the run still draws, the same picture."""
    filename = simulate(str(tmp_path), "vtkhdf")
    with h5py.File(filename, "r") as f:
        assert f["data/dynamic"].shape[0] == 0
    assert picture(filename, 1.0) == picture(siconos_run, 1.0)


def _numbers(text):
    return [float(v) for v in re.split("[ ,]+", text.strip())]


def test_run_and_file_animation_and_coulomb(siconos_run, tmp_path):
    directory = str(tmp_path)
    names = svg.convert(
        siconos_run, animate=True, coulomb=True, times=[0.5], directory=directory
    )
    assert read(names[0]) == read(svg.animation_filename(siconos_run))
    assert read(svg.coulomb_filename(siconos_run, 0.5, directory)) == read(
        svg.coulomb_filename(siconos_run, 0.5)
    )


def test_coulomb_fraction(siconos_run):
    text = read(svg.coulomb_filename(siconos_run, 1.0))
    with h5py.File(siconos_run, "r") as f:
        steps = svg.FileSteps(f["data"])
        contacts = steps.step(steps.resolve(time=1.0)).contacts
    local = contacts[:, svg.CF["reaction_local"]]
    rn, rt, mu = local[:, 0], local[:, 1], contacts[:, svg.CF["mu"]]
    limit = (rn > 0) & (np.abs(rt) >= (1 - svg.SATURATION) * mu * rn)
    assert f"{limit.sum()} of {(rn > 0).sum()} loaded contacts" in text


def test_contact_points_on_outlines(siconos_run):
    """The arrows start on the outline of a disk (or on the slope)."""
    text = read(svg.svg_filename(siconos_run, 1.0))
    circles = [
        _numbers(" ".join(m))
        for m in re.findall(
            r'<circle cx="([-0-9.]+)" cy="([-0-9.]+)" r="([-0-9.]+)"', text
        )
    ]
    starts = re.findall(r'<line x1="([-0-9.]+)" y1="([-0-9.]+)" [^>]*marker-end', text)
    assert starts
    for x, y in starts:
        gaps = [
            abs(np.hypot(float(x) - cx, float(y) - cy) - r) for cx, cy, r in circles
        ]
        assert min(gaps) < 0.05


def test_wall_both_sides(tmp_path, capsys):
    filename = simulate(str(tmp_path), scene="wall")
    capsys.readouterr()
    (name,) = svg.convert(filename, directory=str(tmp_path / "file"))
    text = read(name)
    with h5py.File(filename, "r") as f:
        ids = {name: int(o.attrs["id"]) for name, o in f["data/input"].items()}
    assert f'id="h{-ids["wall"]}_0"' not in text
    assert f'id="h{-ids["floor"]}_0"' in text
    assert f"static {ids['wall']} segment" in capsys.readouterr().err


def _still_contact(text):
    """The cone (e1, p, e2) and the arrow (p, tip) of a picture with one
    contact, in px."""
    cone = re.search(
        r'<path d="M([-0-9.,]+)L([-0-9.,]+)L([-0-9.,]+)" stroke="#e8a33d"', text
    )
    arrow = re.search(
        r'<line x1="([-0-9.]+)" y1="([-0-9.]+)" x2="([-0-9.]+)" y2="([-0-9.]+)" [^>]*marker-end',
        text,
    )
    e1, p, e2 = (np.array(_numbers(g)) for g in cone.groups())
    a = np.array([float(v) for v in arrow.groups()])
    return e1, p, e2, a[:2], a[2:]


def test_cone_saturated_when_sliding(bounce_run, tmp_path):
    """The arrow tip on an edge of the cone at the contacts at the limit,
    inside it at the others; both happen in the run."""
    with h5py.File(bounce_run, "r") as f:
        steps = svg.FileSteps(f["data"])
        found = {True: 0, False: 0}
        for k in range(len(steps.times)):
            contacts = steps.step(k).contacts
            if len(contacts) != 1:
                continue
            (rn, rt, _), mu = contacts[0, svg.CF["reaction_local"]], contacts[0, 0]
            if rn <= 0:
                continue
            saturated = abs(rt) >= (1 - svg.SATURATION) * mu * rn
            (name,) = svg.convert(bounce_run, frames=[k], directory=str(tmp_path))
            e1, p, e2, start, tip = _still_contact(read(name))
            assert np.allclose(start, p, atol=0.01)
            on_edge = min(np.linalg.norm(tip - e1), np.linalg.norm(tip - e2))
            arrow, axis = tip - p, (e1 + e2) / 2 - p
            angle = np.arccos(
                arrow @ axis / np.linalg.norm(arrow) / np.linalg.norm(axis)
            )
            if saturated:
                assert on_edge < 0.05
            else:
                # inside, up to the rounding of the px coordinates
                assert angle < np.arctan(mu) + 5e-4
            found[saturated] += 1
    assert found[True] and found[False]


def test_animation(bounce_run, tmp_path):
    """The animation of the run is the file's; its SMIL lists agree with
    their keyTimes; its frame N is the picture of frame N; the floor has one
    hatch, below it."""
    with h5py.File(bounce_run, "r") as f:
        n = len(svg.FileSteps(f["data"]).times)
    names = svg.convert(bounce_run, animate=True, frames=[7], directory=str(tmp_path))
    animation = read(names[0])
    assert animation == read(svg.animation_filename(bounce_run))
    root = ET.fromstring(animation)
    # the SMIL elements of the frames (the time bar has no keyTimes)
    smil = [
        e
        for e in root.iter()
        if e.tag.endswith(("animate", "animateTransform")) and e.get("keyTimes")
    ]
    assert smil
    for e in smil:
        values = e.get("values").split(";")
        keys = [float(t) for t in e.get("keyTimes").split(";")]
        assert len(values) == len(keys) and keys[0] == 0 and keys == sorted(keys)
    # a value per frame; a linear animation holds the last one to the key
    # time 1
    bodies = next(e for e in root.iter() if e.get("id") == "bodies")
    for e in bodies.iter():
        if e.get("keyTimes") and e.get("attributeName") != "visibility":
            values = e.get("values").split(";")
            if e.get("calcMode") == "linear":
                assert len(values) == n + 1 and values[-1] == values[-2]
                assert e.get("keyTimes").endswith(";1")
            else:
                assert len(values) == n
    body = re.search(
        r'<g id="b1" transform="translate\(([-0-9.,]+)\)"><animateTransform attributeName="transform" type="translate" values="([^"]+)"',
        animation,
    )
    positions = np.array([_numbers(p) for p in body.group(2).split(";")])
    # the contacts of a frame move with their body to its next pose
    followed = 0
    for k in range(n - 1):
        group = next((e for e in root.iter() if e.get("id") == f"contacts{k}"), None)
        if group is None:
            continue
        shift = next(e for e in group.iter() if e.get("type") == "translate")
        end = _numbers(shift.get("values").split(";")[-1])
        assert np.allclose(end, positions[k + 1] - positions[k], atol=0.011)
        followed += 1
    assert followed
    frame = body.group(2).split(";")[7]
    still = read(names[1])
    circle = re.search(
        r'<g id="b1" fill="[^"]*"><circle cx="([-0-9.]+)" cy="([-0-9.]+)"', still
    )
    assert frame == f"{circle.group(1)},{circle.group(2)}"
    assert len(re.findall(r'<path id="h1_0"', animation)) == 1
    stroke = re.search(r'<path id="h1_0" d="M[^l]*l([-0-9.]+),([-0-9.]+)', animation)
    assert float(stroke.group(2)) > 0


def test_animation_controls(bounce_run, tmp_path):
    """--controls adds the buttons and their script, with the start of each
    frame in the animation; without it the animation has no script."""
    with h5py.File(bounce_run, "r") as f:
        n = len(svg.FileSteps(f["data"]).times)
    (plain,) = svg.convert(bounce_run, animate=True, directory=str(tmp_path / "a"))
    (name,) = svg.convert(
        bounce_run, animate=True, controls=True, directory=str(tmp_path / "b")
    )
    assert "<script" not in read(plain)
    text = read(name)
    root = ET.fromstring(text)
    ids = {e.get("id") for e in root.iter()}
    for button in ("first", "previous", "play", "next", "last", "seek", "time"):
        assert f"control-{button}" in ids
    data = re.search(
        r"var ANIMATION = \{starts: \[([^\]]*)\], times: \[([^\]]*)\], "
        r"duration: ([0-9.e+-]+)",
        text,
    )
    starts = _numbers(data.group(1))
    duration = float(data.group(3))
    assert len(starts) == len(_numbers(data.group(2))) == n
    assert starts[0] == 0 and all(np.diff(starts) > 0) and starts[-1] < duration
    # the same animation under the controls
    body = re.compile(r'<g id="bodies".*?\n</g>\n', re.S)
    assert body.search(text).group(0) == body.search(read(plain)).group(0)


def test_restart(tmp_path):
    """A run continued from its file: the converter draws its picture."""
    filename = simulate(str(tmp_path), T=0.5)
    os.remove(svg.svg_filename(filename, 0.5))
    simulate(str(tmp_path), T=1.0, new=False)
    (name,) = svg.convert(filename, times=[1.0], directory=str(tmp_path / "file"))
    assert read(name) == read(svg.svg_filename(filename, 1.0))


# -- the side inference alone -------------------------------------------------


def _contacts(points, normals):
    """cf rows (without time) of contacts with a fixed obstacle."""
    rows = np.zeros((len(points), 25))
    rows[:, svg.CF["point_a"]] = np.pad(points, ((0, 0), (0, 1)))
    rows[:, svg.CF["normal"]] = np.pad(normals, ((0, 0), (0, 1)))
    rows[:, svg.CF["ds1"]] = rows[:, svg.CF["ds2"]] = 1
    return rows


SEGMENT = [((-1, 0), np.array((0.0, 0.0)), np.array((2.0, 0.0)))]


def test_sides_from_contacts():
    sides = svg.Sides(1e-6)
    # bodies above the segment; at the end points, a normal along the
    # segment and the face normal of a polygon pressed on the corner (bench
    # dambreak): no vote
    sides.add_contacts(
        SEGMENT,
        _contacts(
            [(0.5, 0), (1.5, 0), (2, 0), (0, 0)],
            [(0, 1), (0.1, 0.99), (1, 0), (-0.943, -0.332)],
        ),
    )
    (side,) = sides.material(SEGMENT, [], 1.0).values()
    assert np.allclose(side, (0, -1))


def test_sides_both_sides(capsys):
    sides = svg.Sides(1e-6)
    sides.add_contacts(SEGMENT, _contacts([(0.5, 0), (1.5, 0)], [(0, 1), (0, -1)]))
    assert sides.material(SEGMENT, [], 1.0) == {(-1, 0): None}
    assert "bodies on both sides" in capsys.readouterr().err


def test_sides_from_bodies_and_override(capsys):
    sides = svg.Sides(1e-6)
    # no contact: the bodies near the segment, the first set that has some
    (side,) = sides.material(
        SEGMENT, [np.array([[1.0, 5.0]]), np.array([[1.0, -0.5]])], 1.0
    ).values()
    assert np.allclose(side, (0, 1))
    assert sides.material(SEGMENT, [], 1.0) == {(-1, 0): None}
    assert "no contact nor body" in capsys.readouterr().err
    forced = svg.Sides(1e-6, {-1: -1})
    (side,) = forced.material(SEGMENT, [], 1.0).values()
    assert np.allclose(side, (0, -1))


def _box():
    """An open box: left wall down, floor, right wall up; the floor and the
    walls of separate static objects."""
    return [
        ((-1, 0), np.array((0.0, 2.0)), np.array((0.0, 0.0))),
        ((-2, 0), np.array((0.0, 0.0)), np.array((4.0, 0.0))),
        ((-3, 0), np.array((4.0, 0.0)), np.array((4.0, 2.0))),
    ]


def test_sides_chain():
    """Contacts on the floor only: the walls, joined to it, are hatched
    outside."""
    box = _box()
    sides = svg.Sides(1e-6)
    sides.add_contacts(box, _contacts([(2, 0)], [(0, 1)]))
    material = sides.material(box, [], 1.0)
    assert np.allclose(material[(-2, 0)], (0, -1))
    assert np.allclose(material[(-1, 0)], (-1, 0))
    assert np.allclose(material[(-3, 0)], (1, 0))


def test_sides_junction(capsys):
    """A wall standing on the middle of a floor is not in its chain, nor a
    wall at a joint of three segments."""
    floor = [
        ((-1, 0), np.array((0.0, 0.0)), np.array((2.0, 0.0))),
        ((-2, 0), np.array((2.0, 0.0)), np.array((4.0, 0.0))),
        ((-3, 0), np.array((2.0, 0.0)), np.array((2.0, 2.0))),
    ]
    sides = svg.Sides(1e-6)
    sides.add_contacts(floor, _contacts([(1, 0), (3, 0)], [(0, 1), (0, 1)]))
    material = sides.material(floor, [], 1.0)
    assert np.allclose(material[(-1, 0)], (0, -1))
    assert np.allclose(material[(-2, 0)], (0, -1))
    assert material[(-3, 0)] is None
    assert "no contact nor body" in capsys.readouterr().err


def test_parse_times():
    assert svg.parse_times("0.5,1,2e-1, 3.0") == [0.5, 1, 0.2, 3.0]
    assert svg.parse_times("[0.5]") == [0.5]


if __name__ == "__main__":
    run(sys.argv[1], sys.argv[2], sys.argv[3], float(sys.argv[4]), sys.argv[5] == "new")
