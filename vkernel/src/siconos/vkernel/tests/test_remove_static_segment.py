# SPDX-FileCopyrightText: INRIA
#
# SPDX-License-Identifier: Apache-2.0
#
# Siconos is a program dedicated to modeling, simulation and control
# of non smooth dynamical systems.

"""A static segment with a time of death, removed under the vnative backend:
the other segments must still hold the disks.

The removed segment is the last one added, so the first one imported (the
static ids are negative): its points come first, and storage::remove moves
the points of the floor into their slots. Without the z-order sort, the
neighborhood search must still see the moved points where they are, or the
disks fall through the floor.

Run with pytest, nonos and siconos.io in the PYTHONPATH.
"""

import numpy as np
import pytest

import siconos.numerics as sn
from siconos.mechanics.collision.tools import Contactor
from siconos.io.mechanics_run import (
    MechanicsHdf5Runner,
    MechanicsHdf5Runner_run_options,
    RunnerConfig,
)
import nonos

RADIUS = 0.1


def test_remove_static_segment(tmp_path):
    filename = str(tmp_path / "remove_static_segment.hdf5")
    config = RunnerConfig("vnative")

    with MechanicsHdf5Runner(config=config, io_filename=filename) as io:
        io.add_primitive_shape("Disk", "Disk", [RADIUS])
        io.add_Newton_impact_friction_nsl("contact", mu=0.3, e=0)

        # a floor of 10 segments, 20 points each
        for k in range(10):
            io.add_primitive_shape(f"Floor-{k}", "Segment", (k, 0, k + 1, 0))
            io.add_object(f"floor-{k}", [Contactor(f"Floor-{k}")], translation=[0, 0])
        # a wall far from the disks, 200 points: as many as the floor
        io.add_primitive_shape("Wall", "Segment", (11, 0, 11, 10))
        io.add_object(
            "wall", [Contactor("Wall")], translation=[0, 0], time_of_death=0.02
        )

        x = np.linspace(0.25, 9.75, 40)
        io.add_objects(
            "disks",
            [Contactor("Disk")],
            np.column_stack([x, np.full(40, RADIUS)]),
            mass=1.0,
        )

    options = sn.SolverOptions(sn.solver_ids.SICONOS_FRICTION_2D_NSGS)
    vnative_options = nonos.bridge.SpaceFilterOptions()
    vnative_options.neighborhood_radius = 2.5 * RADIUS
    vnative_options.min_radius = 0.5 * RADIUS
    # the sort rebuilds the hash table and would hide a stale one
    vnative_options.sort_frequency = 0

    run_options = MechanicsHdf5Runner_run_options()
    run_options["t0"] = 0
    run_options["T"] = 0.3
    run_options["h"] = 0.005
    run_options["solver_options"] = options
    run_options["vnative_options"] = vnative_options
    run_options["verbose"] = False

    with MechanicsHdf5Runner(mode="r+", config=config, io_filename=filename) as io:
        io.run(run_options)
        y = np.array([io._nsds.dynamicalSystem(i)._bdy.q()[1] for i in range(1, 41)])

    # resting: 0.3 s of free fall would be 0.44 below
    assert y.min() > 0.5 * RADIUS, y.min()


@pytest.mark.parametrize("box_added_first", [True, False])
def test_remove_static_box2d(tmp_path, box_added_first):
    """A Box2d is four segments removed one after the other. The Python
    handles hold their index by value: when the four segments are the last
    ones (box added first, so imported last), a removal moves one of them
    into the freed slot and its handle is stale.

    The disks rest on the top of the box; once the box is removed they must
    fall onto the floor below, not on a leftover side of the box."""
    filename = str(tmp_path / "remove_static_box2d.hdf5")
    config = RunnerConfig("vnative")
    height = 0.5

    with MechanicsHdf5Runner(config=config, io_filename=filename) as io:
        io.add_primitive_shape("Disk", "Disk", [RADIUS])
        io.add_Newton_impact_friction_nsl("contact", mu=0.3, e=0)

        def add_box():
            # a flat box: 10.5 wide, 0.2 high, its top at height
            io.add_primitive_shape(
                "Box", "Box2d", (10.5, 0.2), insideMargin=0.0, outsideMargin=0.0
            )
            io.add_object(
                "box",
                [Contactor("Box")],
                translation=[5, height - 0.1],
                time_of_death=0.02,
            )

        if box_added_first:
            add_box()
        for k in range(10):
            io.add_primitive_shape(f"Floor-{k}", "Segment", (k, 0, k + 1, 0))
            io.add_object(f"floor-{k}", [Contactor(f"Floor-{k}")], translation=[0, 0])
        if not box_added_first:
            add_box()

        x = np.linspace(0.25, 9.75, 40)
        io.add_objects(
            "disks",
            [Contactor("Disk")],
            np.column_stack([x, np.full(40, height + RADIUS)]),
            mass=1.0,
        )

    options = sn.SolverOptions(sn.solver_ids.SICONOS_FRICTION_2D_NSGS)
    vnative_options = nonos.bridge.SpaceFilterOptions()
    vnative_options.neighborhood_radius = 2.5 * RADIUS
    vnative_options.min_radius = 0.5 * RADIUS
    vnative_options.sort_frequency = 0

    run_options = MechanicsHdf5Runner_run_options()
    run_options["t0"] = 0
    # 0.5 of free fall takes 0.32 s
    run_options["T"] = 0.6
    run_options["h"] = 0.005
    run_options["solver_options"] = options
    run_options["vnative_options"] = vnative_options
    run_options["verbose"] = False

    with MechanicsHdf5Runner(mode="r+", config=config, io_filename=filename) as io:
        io.run(run_options)
        y = np.array([io._nsds.dynamicalSystem(i)._bdy.q()[1] for i in range(1, 41)])

    # on the floor: y = RADIUS; on the bottom of the box: height - 0.2 + RADIUS
    assert y.min() > 0.5 * RADIUS, y.min()
    assert y.max() < 2 * RADIUS, y.max()
