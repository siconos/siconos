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

has_mechanics_run = False

try:
    import siconos.io.mechanics_run as smrun
    import siconos.mechanics.collision.tools as smct

    has_mechanics_run = True
except Exception:
    pass


# @pytest.mark.skipif(not has_mechanics_run, reason="Only if mechanics is available")
def test_create_h5run():
    outputfile_name = "mytest.h5"
    # Just test creation and proper writing for the output file
    with smrun.MechanicsHdf5Runner(io_filename=outputfile_name) as myio:
        myio.print_verbose("Build mechanics_hdf5")
        # Definition of a cube as a convex shape
        myio.add_convex_shape(
            "Cube",
            [
                (-1.0, 1.0, -1.0),
                (-1.0, -1.0, -1.0),
                (-1.0, -1.0, 1.0),
                (-1.0, 1.0, 1.0),
                (1.0, 1.0, 1.0),
                (1.0, 1.0, -1.0),
                (1.0, -1.0, -1.0),
                (1.0, -1.0, 1.0),
            ],
        )

        # Alternative to the previous convex shape definition.
        # io.add_primitive_shape('Cube1', 'Box', (2, 2, 2))

        # Definition of the ground shape
        myio.add_primitive_shape("Ground", "Box", (100, 100, 0.5))

        # Definition of a non smooth law. As no group ids are specified it
        # is between contactors of group id 0.
        myio.add_Newton_impact_friction_nsl("contact", mu=0.3)

        # The cube object made with an unique Contactor : the cube shape.
        # As a mass is given, it is a dynamic system involved in contact
        # detection and in the simulation.  With no group id specified the
        # Contactor belongs to group 0
        myio.add_object(
            "cube",
            [smct.Contactor("Cube")],
            translation=[0, 0, 2],
            velocity=[10, 0, 0, 1, 1, 1],
            mass=1,
        )

        # the ground object made with the ground shape. As the mass is
        # not given, it is a static object only involved in contact
        # detection.
        myio.add_object(
            "ground",
            [smct.Contactor("Ground")],
            translation=[0, 0, 0],
        )

    fpath = Path(outputfile_name)
    assert fpath.is_file(), f"Error : file '{outputfile_name}' does not exist!"

    # Check if the file is a valid HDF5 file
    try:
        with h5py.File(outputfile_name, "r") as f:
            # Verify that the 'data' group exists
            assert (
                "data" in f
            ), "Error: The HDF5 file does not contain the 'data' group."

            # Verify the 'dimension' attribute exists
            assert (
                "dimension" in f.attrs
            ), "Error: The 'dimension' attribute is missing."

            # Check results of add_convex_shape
            # Verify if the 'dimension' attribute exists
            assert (
                "dimension" in f.attrs
            ), "Error: The 'dimension' attribute is missing."

            assert f.attrs["dimension"] == 3
            assert "Cube" in f["data"]["ref"]
            cube = f["data"]["ref"]["Cube"]
            assert cube.shape == (8, 3)

            # Then add_primitive_shape
            assert "Ground" in f["data"]["ref"]
            ground = f["data"]["ref"]["Ground"]
            assert ground.shape == (1, 3)
            assert ground.attrs["primitive"] == "Box"

            # add_Newton_impact_friction_nsl
            assert "contact" in f["data"]["nslaws"]
            contact = f["data"]["nslaws"]["contact"]
            assert contact.attrs["type"] == "NewtonImpactFrictionNSL"
            assert contact.attrs["mu"] == 0.3

            # add_object
            assert "cube" in f["data"]["input"]
            assert "ground" in f["data"]["input"]
            cube_obj = f["data"]["input"]["cube"]
            assert (cube_obj.attrs["translation"] == [0, 0, 2]).all()
            assert (cube_obj.attrs["velocity"] == [10, 0, 0, 1, 1, 1]).all()
            assert "Cube-0" in cube_obj

    except OSError:
        pytest.fail(f"Error: '{outputfile_name}' is not a valid HDF5 file!")

    # Optional: Remove the file after the test to keep the workspace clean
    # os.remove(outputfile_name)
