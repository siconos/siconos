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
"""Tests of siconos_filter on the output of a simulation: a pendulum (a ball
on a pivot joint) that hits a wall.

The filtered files are read back with format_hdf5 (as vview does): every table
must keep the rows of the kept steps, and its step index must stay consistent.
"""

import shutil
import subprocess
import sys
from pathlib import Path

import h5py
import numpy as np
import pytest

has_mechanics_run = False
try:
    import siconos.io
    import siconos.io.mechanics_run as smrun
    import siconos.mechanics
    import siconos.mechanics.collision.tools as smct
    from siconos.io import format_hdf5

    has_mechanics_run = True
except ImportError:
    pass

has_bullet = has_mechanics_run and siconos.mechanics.have_bullet


def filter_script():
    """siconos_filter next to the siconos.io package (build) or in the PATH (install)."""
    script = Path(siconos.io.__file__).parent / "siconos_filter"
    if script.is_file():
        return str(script)
    return shutil.which("siconos_filter")


pytestmark = [
    pytest.mark.skipif(
        not has_bullet, reason="Only if mechanics with bullet is available"
    ),
    pytest.mark.skipif(
        has_bullet and filter_script() is None,
        reason="Only if siconos_filter is available",
    ),
]

T, h = 1.5, 5e-3
OBJECTS_TABLES = ("dynamic", "velocities", "static")  # ds_id in column 1


@pytest.fixture(scope="module")
def pendulum(tmp_path_factory):
    """Output of the simulation: the ball, length 1, released horizontally,
    swings down and hits the wall on the other side."""
    filename = str(tmp_path_factory.mktemp("filter") / "pendulum.hdf5")
    config = smrun.RunnerConfig(backend="bullet")
    with smrun.MechanicsHdf5Runner(io_filename=filename, config=config) as io:
        io.add_primitive_shape("Ball", "Sphere", (0.1,))
        io.add_primitive_shape("Wall", "Box", (0.1, 1.0, 1.0))
        io.add_Newton_impact_friction_nsl("contact", mu=0.3, e=0.5)
        io.add_object("ball", [smct.Contactor("Ball")], translation=[1, 0, 0], mass=1)
        # contact with the ball at x = -0.5 (its center), z = -0.87
        io.add_object("wall", [smct.Contactor("Wall")], translation=[-0.65, 0, -1])
        io.add_joint(
            "pivot",
            "ball",
            None,
            [[0, 0, 0]],
            [[0, 1, 0]],
            "PivotJointR",
            absolute=True,
        )

    run_options = smrun.MechanicsHdf5Runner_run_options()
    run_options["T"] = T
    run_options["h"] = h
    with smrun.MechanicsHdf5Runner(
        io_filename=filename, mode="r+", config=config
    ) as io:
        io.run(run_options)

    with h5py.File(filename, "r") as f:
        assert f["data/cf"].shape[0] > 0, "the pendulum does not hit the wall"
    return filename


def run_filter(fin, fout, *options):
    return subprocess.run(
        [sys.executable, filter_script(), str(fin), str(fout), *options],
        capture_output=True,
        text=True,
        check=False,  # the return code is checked by the tests
    )


def filtered(fin, tmp_path, *options):
    fout = tmp_path / "out.hdf5"
    result = run_filter(fin, fout, *options)
    assert result.returncode == 0, result.stdout + result.stderr
    return fout


def object_ids(name, rows):
    """Columns of the objects of the rows of a table, None if it has none."""
    if name == "cf":
        return rows[:, -2:]  # ds1, ds2
    if name == "cf_info":
        return rows[:, 2:4]  # ds1, ds2
    if name in OBJECTS_TABLES:
        return rows[:, 1:2]
    return None


def check(fin, fout, kept_times=None, excluded_ids=(), single=False):
    """The output is read with the format, its steps are the kept steps of the
    input (all of them by default), with the rows of the input that do not
    reference an excluded object; the step index is consistent with its table."""
    fa_in = format_hdf5.FormatAdaptor(str(fin), "r")
    fa_out = format_hdf5.FormatAdaptor(str(fout), "r")
    try:
        fa_out.open()
        assert dict(fa_in.hdf5.attrs) == dict(fa_out.hdf5.attrs)
        indexed = "step_index" in fa_in.data
        for name, spec in fa_in.data_format.items():
            if "<" in name or "dtype" not in spec or name not in fa_in.data:
                continue
            steps_in, steps_out = fa_in.rows(name), fa_out.rows(name)
            if steps_in is None:
                continue
            times_in, starts_in, ends_in = steps_in
            times_out, starts_out, ends_out = steps_out
            table_in, table_out = fa_in.data[name][:], fa_out.data[name][:]
            if indexed:
                index = fa_out.data["step_index"][name]
                assert index.dtype == np.float64, name
                assert np.array_equal(starts_out[1:], ends_out[:-1]), name
                assert (ends_out[-1] if len(ends_out) else 0) == table_out.shape[0], (
                    name
                )
            if single:
                assert table_out.dtype == np.float32, name

            keep = np.ones(len(times_in), dtype=bool)
            if kept_times is not None:
                keep = np.isin(
                    times_in.astype("f4"), np.asarray(kept_times, dtype="f4")
                )
            ids = object_ids(name, table_in)
            left = np.ones(table_in.shape[0], dtype=bool)
            if ids is not None and len(excluded_ids):
                left = ~np.isin(ids, excluded_ids).any(axis=1)
            if not indexed:  # steps from column t: a step left empty disappears
                keep &= [left[s:e].any() for s, e in zip(starts_in, ends_in)]

            assert np.array_equal(
                times_out.astype("f4"), times_in[keep].astype("f4")
            ), name
            for s_in, e_in, s_out, e_out in zip(
                starts_in[keep], ends_in[keep], starts_out, ends_out
            ):
                rows = table_in[s_in:e_in][left[s_in:e_in]]
                assert np.allclose(rows, table_out[s_out:e_out], equal_nan=True), name
    finally:
        fa_in.close()
        fa_out.close()


def test_copy(pendulum, tmp_path):
    fout = filtered(pendulum, tmp_path)
    check(pendulum, fout)
    with h5py.File(pendulum, "r") as fin, h5py.File(fout, "r") as f:
        assert f.attrs["version_major"] == 1
        assert f["format"][()] == fin["format"][()]
        assert "pivot" in f["data/joints"]


def test_single_gzip(pendulum, tmp_path):
    fout = filtered(pendulum, tmp_path, "--single", "--gzip")
    check(pendulum, fout, single=True)


@pytest.mark.parametrize("start, end", [(0.5, 1.0), (None, 0.0), (100.0, None)])
def test_time_window(pendulum, tmp_path, start, end):
    options = []
    with h5py.File(pendulum, "r") as f:
        times = np.unique(f["data/dynamic"][:, 0])
    kept = times
    if start is not None:
        options += ["--start", str(start)]
        kept = kept[kept.astype("f8") >= start]
    if end is not None:
        options += ["--end", str(end)]
        kept = kept[kept.astype("f8") <= end]
    fout = filtered(pendulum, tmp_path, *options)
    check(pendulum, fout, kept_times=kept)


def test_interval(pendulum, tmp_path):
    interval = 0.05
    fout = filtered(pendulum, tmp_path, "--interval", str(interval))
    with h5py.File(fout, "r") as f:
        kept = np.unique(f["data/dynamic"][:, 0])
    assert 0 < len(kept) < T / h
    assert np.all(np.diff(kept) > interval / 2)
    check(pendulum, fout, kept_times=kept)


def test_exclude(pendulum, tmp_path):
    with h5py.File(pendulum, "r") as f:
        ball = f["data/input/ball"].attrs["id"]
    fout = filtered(pendulum, tmp_path, "--exclude", "ball")
    check(pendulum, fout, excluded_ids=[ball])
    with h5py.File(fout, "r") as f:
        assert "ball" not in f["data/input"]
        assert "pivot" not in f["data/joints"]  # names the ball
        assert f["data/cf"].shape[0] == 0


def test_attr(pendulum, tmp_path):
    fout = filtered(pendulum, tmp_path, "--attr", "nslaws.contact.mu=0.7")
    check(pendulum, fout)
    with h5py.File(fout, "r") as f:
        assert f["data/nslaws/contact"].attrs["mu"] == 0.7


def test_file_without_version(pendulum, tmp_path):
    """A file without version attributes (nor format, step index) is a v1.0 file."""
    old = tmp_path / "old.hdf5"
    shutil.copy(pendulum, old)
    with h5py.File(old, "r+") as f:
        for name in ("version_major", "version_minor"):
            del f.attrs[name]
        del f["format"]
        del f["data/step_index"]
    fout = filtered(old, tmp_path, "--start", "0.5", "--exclude", "wall")
    with h5py.File(old, "r") as f:
        wall = f["data/input/wall"].attrs["id"]
    with h5py.File(fout, "r") as f:
        assert "version_major" not in f.attrs
        times = np.unique(f["data/dynamic"][:, 0])
        assert np.all(times.astype("f8") >= 0.5)
        assert not np.isin(f["data/static"][:, 1], [wall]).any()


def test_other_version_refused(pendulum, tmp_path):
    v2 = tmp_path / "v2.hdf5"
    shutil.copy(pendulum, v2)
    with h5py.File(v2, "r+") as f:
        f.attrs["version_major"] = 2
    result = run_filter(v2, tmp_path / "out.hdf5")
    assert result.returncode != 0
    assert not (tmp_path / "out.hdf5").exists()
