"""
Standalone pytest for FEM strain/stress tensor HDF5 output.

This test exercises the vkernel FEM path end-to-end:
  1. Run ``fem.py`` on a small Gmsh mesh.
  2. Open the produced ``fem.hdf5``.
  3. Verify that strain and stress datasets were written with the
     expected names, attributes, and layout.

It is intentionally kept outside CMake/CTest so it can be run directly
with ``pytest`` once the Siconos Python environment is available.
"""

from __future__ import annotations

import os
import subprocess
import sys

import h5py
import numpy as np
import pytest

# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

# The test files live next to this one.
_HERE = os.path.dirname(os.path.abspath(__file__))

# Absolute path to the tiny Gmsh mesh shipped with the vkernel tests.
MESH_PATH = os.path.join(_HERE, "mesh_N5.msh")

# ``fem.py`` lives next to this test file.
FEM_SCRIPT = os.path.join(_HERE, "fem.py")

# Analytical affine-field regression check (see module docstring there).
ANALYTICAL_SCRIPT = os.path.join(_HERE, "fem_tensor_analytical.py")

# Siconos Python bindings are not installed system-wide: they are built
# in-tree and exposed via PYTHONPATH and LD_LIBRARY_PATH, from the build
# directory given by SICONOS_BUILD_DIR.
_BUILD_DIR = os.environ.get("SICONOS_BUILD_DIR", "")

_SICONOS_PYTHONPATHS = [
    os.path.join(_HERE, "..", "py"),
    os.path.join(_BUILD_DIR, "python"),
]

# The shared libraries produced by the build are not on the system loader
# path by default.
_SICONOS_LIBRARY_PATHS = [
    os.path.join(_BUILD_DIR, d)
    for d in (
        "numerics",
        "mechanics",
        "kernel",
        "io",
        "externals",
        "vkernel/soa",
        "control",
        "vkernel/soa/extern/CompactNSearch",
    )
]


def _fem_subprocess_env():
    """Environment (PYTHONPATH/LD_LIBRARY_PATH) needed to run the in-tree
    Siconos Python bindings of SICONOS_BUILD_DIR."""
    if not _BUILD_DIR:
        pytest.skip("SICONOS_BUILD_DIR is not set")
    env = dict(os.environ)
    existing_pythonpath = env.get("PYTHONPATH", "")
    pythonpath_extra = ":".join(_SICONOS_PYTHONPATHS)
    env["PYTHONPATH"] = (
        (pythonpath_extra + ":" + existing_pythonpath)
        if existing_pythonpath
        else pythonpath_extra
    )

    existing_ld = env.get("LD_LIBRARY_PATH", "")
    ld_extra = ":".join(_SICONOS_LIBRARY_PATHS)
    env["LD_LIBRARY_PATH"] = (ld_extra + ":" + existing_ld) if existing_ld else ld_extra
    return env


# The Siconos Python extensions must be loaded by the interpreter they were
# compiled for: SICONOS_PYTHON, or the interpreter running the tests.
_PYTHON = os.environ.get("SICONOS_PYTHON", sys.executable)


def _run_fem_simulation(tmp_path, mesh_path: str) -> None:
    """Run ``fem.py`` in an isolated temp directory.

    The script writes ``fem.hdf5`` in its current working directory, so we
    cd into ``tmp_path`` beforehand to avoid polluting the source tree.
    """
    # Run in tmp_path so fem.hdf5 is created there.
    subprocess.run(
        [_PYTHON, FEM_SCRIPT, mesh_path],
        cwd=tmp_path,
        env=_fem_subprocess_env(),
        check=True,
    )


def _run_analytical_check(tmp_path, mesh_path: str) -> subprocess.CompletedProcess:
    """Run ``fem_tensor_analytical.py`` in an isolated temp directory.

    Unlike ``_run_fem_simulation``, this does not assert success directly
    (``check=False``) so the caller can print the script's diagnostic
    output on failure.
    """
    return subprocess.run(
        [_PYTHON, ANALYTICAL_SCRIPT, mesh_path],
        cwd=tmp_path,
        env=_fem_subprocess_env(),
        check=False,
        capture_output=True,
        text=True,
    )


# ---------------------------------------------------------------------------
# Tests
# ---------------------------------------------------------------------------


def test_fem_tensor_datasets_exist(tmp_path):
    """The HDF5 file must contain both strain and stress datasets.

    Dataset naming convention: ``fem_<tensor>_<ds_id>`` where ``ds_id``
    is the dynamical-system id assigned by the mechanics runner.
    """
    _run_fem_simulation(tmp_path, MESH_PATH)

    hdf5_path = tmp_path / "fem.hdf5"
    assert hdf5_path.is_file(), "fem.py did not produce fem.hdf5"

    with h5py.File(hdf5_path, "r") as f:
        # The mechanics runner stores per-DS outputs under the ``data`` group.
        assert "data" in f, "HDF5 file missing 'data' group"

        # At least one FEM ds_id should have been registered.
        epsilon_groups = [k for k in f["data"].keys() if k.startswith("fem_epsilon_")]
        sigma_groups = [k for k in f["data"].keys() if k.startswith("fem_sigma_")]

        assert epsilon_groups, "No fem_epsilon_* datasets found"
        assert sigma_groups, "No fem_sigma_* datasets found"

        # Strain and stress should be written for the same ds_id(s).
        epsilon_ids = {k.split("_")[-1] for k in epsilon_groups}
        sigma_ids = {k.split("_")[-1] for k in sigma_groups}
        assert epsilon_ids == sigma_ids, (
            f"Strain/stress ds_ids differ: epsilon={epsilon_ids}, sigma={sigma_ids}"
        )


def test_fem_tensor_dataset_attributes(tmp_path):
    """Each tensor dataset must carry the expected attributes.

    ``output_fem_epsilon`` / ``output_fem_sigma`` in ``mechanics_run.py``
    set:
      - ``tensor_type = "symmetric_2d"``
      - ``components``  (3 strings: exx/eyy/exy or sxx/syy/sxy)
      - ``num_elements`` (number of T3 elements in the mesh)
    """
    _run_fem_simulation(tmp_path, MESH_PATH)

    with h5py.File(tmp_path / "fem.hdf5", "r") as f:
        data = f["data"]

        for ds_name in data.keys():
            if not ds_name.startswith(("fem_epsilon_", "fem_sigma_")):
                continue

            ds = data[ds_name]

            # Mandatory scalar attributes written by mechanics_run.py.
            assert "tensor_type" in ds.attrs, f"{ds_name} missing tensor_type"
            assert ds.attrs["tensor_type"] == "symmetric_2d", (
                f"{ds_name} has unexpected tensor_type={ds.attrs['tensor_type']}"
            )

            assert "components" in ds.attrs, f"{ds_name} missing components"
            components = list(ds.attrs["components"])
            assert len(components) == 3, f"{ds_name} components length != 3"

            assert "num_elements" in ds.attrs, f"{ds_name} missing num_elements"
            num_elements = int(ds.attrs["num_elements"])
            assert num_elements > 0, f"{ds_name} has non-positive num_elements"


def test_fem_tensor_dataset_layout(tmp_path):
    """Dataset columns must follow the flat per-element convention.

    mechanics_run.py writes rows as::

        [time, ds_id, exx_0, eyy_0, exy_0, exx_1, eyy_1, exy_1, ...]

    so ``shape[1] == 2 + 3 * num_elements`` and every row has the same
    ``ds_id`` in column 1.
    """
    _run_fem_simulation(tmp_path, MESH_PATH)

    with h5py.File(tmp_path / "fem.hdf5", "r") as f:
        data = f["data"]

        for ds_name in data.keys():
            if not ds_name.startswith(("fem_epsilon_", "fem_sigma_")):
                continue

            ds = data[ds_name]
            num_elements = int(ds.attrs["num_elements"])
            expected_cols = 2 + 3 * num_elements

            # At least one output row must exist.
            assert ds.shape[0] > 0, f"{ds_name} has no time rows"
            assert ds.shape[1] == expected_cols, (
                f"{ds_name} columns={ds.shape[1]}, expected {expected_cols}"
            )

            # Every row must carry the same ds_id in column 1.
            ds_ids = ds[:, 1].astype(int)
            expected_ds_id = int(ds_name.split("_")[-1])
            assert np.all(ds_ids == expected_ds_id), (
                f"{ds_name} has inconsistent ds_id values"
            )


def test_fem_tensor_values_are_finite(tmp_path):
    """Tensor entries should be finite numbers at every output step.

    A non-finite value here indicates a crash inside the FEM strain/stress
    computation or a broken write path.
    """
    _run_fem_simulation(tmp_path, MESH_PATH)

    with h5py.File(tmp_path / "fem.hdf5", "r") as f:
        data = f["data"]

        for ds_name in data.keys():
            if not ds_name.startswith(("fem_epsilon_", "fem_sigma_")):
                continue

            values = data[ds_name][:, 2:]  # skip time + ds_id columns
            assert np.all(np.isfinite(values)), (
                f"{ds_name} contains non-finite tensor values"
            )


def test_fem_tensor_matches_analytical_affine_field(tmp_path):
    """Strain/stress must be numerically correct, not just finite.

    Imposes a synthetic affine displacement field on the mesh (see
    ``fem_tensor_analytical.py``) and checks that every T3 element reports
    the exact same, analytically-known strain/stress. A CST element is
    exact for any affine field regardless of element size, so any
    element-size-dependent scaling bug in the strain-displacement matrix
    (e.g. reusing an area-scaled assembly matrix instead of the plain
    kinematic one) would show up here as per-element mismatches, even
    though the earlier finiteness/shape checks above would still pass.
    """
    result = _run_analytical_check(tmp_path, MESH_PATH)
    assert result.returncode == 0, (
        f"analytical tensor check failed:\n{result.stdout}\n{result.stderr}"
    )
