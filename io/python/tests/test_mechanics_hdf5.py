# Siconos is a program dedicated to modeling, simulation and control
# of non smooth dynamical systems.
#
# SPDX-FileCopyrightText: INRIA
# SPDX-License-Identifier: Apache-2.0
#
import os
from pathlib import Path

import h5py
import pytest

has_mechanics_io = False
try:
    import siconos.io.mechanics_hdf5 as sh5

    has_mechanics_io = True
except ImportError:
    pass


@pytest.mark.skipif(not has_mechanics_io, reason="Only if mechanics is available")
def test_create_h5():
    outputfile_name = "mytest.h5"
    # Just test creation and proper writing for the output file
    with sh5.MechanicsHdf5(io_filename=outputfile_name) as myio:
        myio.print_verbose("Build mechanics_hdf5")

    fpath = Path(outputfile_name)
    assert fpath.is_file(), f"Error : file '{outputfile_name}' does not exist!"

    # Check if the file is a valid HDF5 file
    try:
        with h5py.File(outputfile_name, "r") as f:
            # Verify that the 'data' group exists
            assert "data" in f, (
                "Error: The HDF5 file does not contain the 'data' group."
            )

            # Verify that 'static' and 'velocities' datasets exist
            assert "static" in f["data"], "The 'static' dataset is missing."
            assert "velocities" in f["data"], (
                "Error: The 'velocities' dataset is missing."
            )

            # Check the dimensions of the datasets
            static_shape = f["data"]["static"].shape
            velocities_shape = f["data"]["velocities"].shape
            assert static_shape == (
                0,
                9,
            ), "Error: 'static' does not have the correct (n x m) shape."
            assert velocities_shape == (
                0,
                8,
            ), "Error: 'velocities' does not have the correct (n x m) shape."

    except OSError:
        pytest.fail(f"Error: '{outputfile_name}' is not a valid HDF5 file!")

    os.remove(outputfile_name)


if __name__ == "__main__":
    test_create_h5()
