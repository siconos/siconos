"""add_objects with per-grain radii and masses, imported with the vnative
backend.

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

SHAPE_RADIUS = 0.1
NEIGHBORHOOD_RADIUS = 2.1 * SHAPE_RADIUS
GAP = 0.02


def grid(n, x0, spacing=1.0):
    return np.column_stack([x0 + spacing * np.arange(n), np.zeros(n)])


def test_add_objects_polydisperse(tmp_path):
    filename = str(tmp_path / "add_objects.hdf5")
    config = RunnerConfig("vnative")

    # the last radius is bigger than NEIGHBORHOOD_RADIUS / 2: several
    # points on its circle for the neighborhood search
    radii = np.array([0.05, 0.1, 0.2, 0.5])
    masses = np.array([1.0, 2.0, 3.0, 4.0])
    density = 10.0

    with MechanicsHdf5Runner(config=config, io_filename=filename) as io:
        io.add_primitive_shape("Disk", "Disk", [SHAPE_RADIUS])
        io.add_Newton_impact_friction_nsl("contact", mu=0.3, e=0)

        # ids 1-4: per-grain radii, mass from the density. The grains
        # touch each other, and the big one comes at 1 m/s from 0.02 away:
        # its contact is only seen through the points on its circle.
        x = np.concatenate([[0.0], np.cumsum(radii[:-1] + radii[1:])])
        x[-1] += GAP
        velocities = np.zeros((4, 3))
        velocities[-1, 0] = -1.0
        io.add_objects(
            "by_density",
            [Contactor("Disk")],
            np.column_stack([x, np.zeros(4)]),
            velocities=velocities,
            radii=radii,
            density=density,
        )
        # ids 5-8: shared radius different from the shape one, shared mass
        io.add_objects(
            "scalar", [Contactor("Disk")], grid(4, 10.0), radii=0.3, mass=1.5
        )
        # ids 9-12: per-grain radii and masses
        io.add_objects(
            "per_grain", [Contactor("Disk")], grid(4, 20.0), radii=radii, mass=masses
        )
        # ids 13-16: radius of the shape
        io.add_objects("default", [Contactor("Disk")], grid(4, 30.0), mass=1.0)

        with pytest.raises(ValueError):
            io.add_objects(
                "both", [Contactor("Disk")], grid(4, 40.0), mass=1.0, density=1.0
            )
        with pytest.raises(ValueError):
            io.add_objects(
                "bad_radii", [Contactor("Disk")], grid(4, 50.0), radii=radii[:3]
            )

    options = sn.SolverOptions(sn.solver_ids.SICONOS_FRICTION_2D_NSGS)
    vnative_options = nonos.bridge.SpaceFilterOptions()
    vnative_options.neighborhood_radius = NEIGHBORHOOD_RADIUS
    vnative_options.min_radius = 0.5 * SHAPE_RADIUS

    run_options = MechanicsHdf5Runner_run_options()
    run_options["t0"] = 0
    run_options["T"] = 0.05
    run_options["h"] = 0.005
    run_options["solver_options"] = options
    run_options["vnative_options"] = vnative_options
    run_options["verbose"] = False

    with MechanicsHdf5Runner(mode="r+", config=config, io_filename=filename) as io:
        io.run(run_options)

        def check(ids, expected_radii, expected_masses):
            for i, r, m in zip(ids, expected_radii, expected_masses):
                body = io._nsds.dynamicalSystem(i)._bdy
                assert body.shape().radius() == pytest.approx(r)
                assert list(body.mass_matrix()) == pytest.approx([m, m, m * r * r / 2])

        check(range(1, 5), radii, density * np.pi * radii**2)

        # 0.05 at 1 m/s: a missed contact would end 0.03 inside. A
        # seen one may still end up to one step inside, h * 1 m/s.
        q2 = io._nsds.dynamicalSystem(3)._bdy.q()
        q3 = io._nsds.dynamicalSystem(4)._bdy.q()
        assert q3[0] - q2[0] - (radii[2] + radii[3]) > -run_options["h"]
        check(range(5, 9), [0.3] * 4, [1.5] * 4)
        check(range(9, 13), radii, masses)
        check(range(13, 17), [SHAPE_RADIUS] * 4, [1.0] * 4)
