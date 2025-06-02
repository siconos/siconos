# Siconos is a program dedicated to modeling, simulation and control
# of non smooth dynamical systems.
#
# Copyright 2025 INRIA.
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
from pathlib import Path
import h5py
import pytest
import os

has_mechanics_io = False
try:
    import siconos.io.mechanics_hdf5 as sh5

    has_mechanics_io = True
except Exception:
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
            assert (
                "data" in f
            ), "Error: The HDF5 file does not contain the 'data' group."

            # Verify that 'static' and 'velocities' datasets exist
            assert "static" in f["data"], "The 'static' dataset is missing."
            assert (
                "velocities" in f["data"]
            ), "Error: The 'velocities' dataset is missing."

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
