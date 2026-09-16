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

import subprocess
import sys

import h5py
import numpy as np
import pytest

# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

# Absolute path to the tiny Gmsh mesh shipped with the vkernel tests.
MESH_PATH = (
    "/home/maurice/wkt/siconos/main-devel-constraint-tensor/"
    "siconos/vkernel/soa/src/siconos/tests/mesh_N5.msh"
)

# ``fem.py`` lives next to this test file.
FEM_SCRIPT = (
    "/home/maurice/wkt/siconos/main-devel-constraint-tensor/"
    "siconos/vkernel/soa/src/siconos/tests/fem.py"
)

# Siconos Python bindings are not installed system-wide in this workspace.
# They are built in-tree and exposed via PYTHONPATH, matching the workflow
# documented in STATE.org / GUIX-INSTALL.org.
_SICONOS_PYTHONPATHS = [
    "/home/maurice/wkt/siconos/main-devel-constraint-tensor/"
    "siconos/vkernel/soa/src/siconos/py",
    "/scratch/maurice/agent-atlas/build/Debug-constraint-tensor-guix/python",
]

# The shared libraries produced by the build are not on the system loader
# path by default. STATE.org prescribes the exact LD_LIBRARY_PATH needed
# to run the Python bindings in-tree.
_SICONOS_LIBRARY_PATHS = [
    "/scratch/maurice/agent-atlas/build/Debug-constraint-tensor-guix/numerics",
    "/scratch/maurice/agent-atlas/build/Debug-constraint-tensor-guix/mechanics",
    "/scratch/maurice/agent-atlas/build/Debug-constraint-tensor-guix/kernel",
    "/scratch/maurice/agent-atlas/build/Debug-constraint-tensor-guix/io",
    "/scratch/maurice/agent-atlas/build/Debug-constraint-tensor-guix/externals",
    "/scratch/maurice/agent-atlas/build/Debug-constraint-tensor-guix/vkernel/soa",
    "/scratch/maurice/agent-atlas/build/Debug-constraint-tensor-guix/control",
    "/scratch/maurice/agent-atlas/build/Debug-constraint-tensor-guix/vkernel/soa/extern/CompactNSearch",
]


def _run_fem_simulation(tmp_path, mesh_path: str) -> None:
    """Run ``fem.py`` in an isolated temp directory.

    The script writes ``fem.hdf5`` in its current working directory, so we
    cd into ``tmp_path`` beforehand to avoid polluting the source tree.
    """
    # The Siconos Python extensions in this build tree are compiled for
    # Python 3.11, so we must use that interpreter explicitly rather than
    # whatever ``sys.executable`` points to.
    _PYTHON = (
        "/home/maurice/wkt/siconos/main-devel-constraint-tensor/"
        "siconos/.venv311/bin/python"
    )

    env = dict(__import__("os").environ)
    existing_pythonpath = env.get("PYTHONPATH", "")
    pythonpath_extra = ":".join(_SICONOS_PYTHONPATHS)
    env["PYTHONPATH"] = (
        (pythonpath_extra + ":" + existing_pythonpath) if existing_pythonpath else pythonpath_extra
    )

    existing_ld = env.get("LD_LIBRARY_PATH", "")
    ld_extra = ":".join(_SICONOS_LIBRARY_PATHS)
    env["LD_LIBRARY_PATH"] = (
        (ld_extra + ":" + existing_ld) if existing_ld else ld_extra
    )

    # Run in tmp_path so fem.hdf5 is created there.
    subprocess.run(
        [_PYTHON, FEM_SCRIPT, mesh_path],
        cwd=tmp_path,
        env=env,
        check=True,
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
