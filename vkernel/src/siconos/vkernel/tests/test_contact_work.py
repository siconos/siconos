# SPDX-FileCopyrightText: 2026 INRIA
#
# SPDX-License-Identifier: Apache-2.0
#
# Siconos is a program dedicated to modeling, simulation and control
# of non smooth dynamical systems.

"""cf_work, the contact work table, with the vnative backend.

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


def test_contact_work(tmp_path):
    filename = str(tmp_path / "contact_work.hdf5")
    config = RunnerConfig("vnative")

    # three touching grains, a fourth 0.02 away coming at 1 m/s: a
    # neighbour without contact (an inactive interaction) for the first
    # steps, then an impact
    n = 4
    x = 2 * RADIUS * np.arange(n)
    x[-1] += 0.02
    velocities = np.zeros((n, 3))
    velocities[-1, 0] = -1.0

    with MechanicsHdf5Runner(config=config, io_filename=filename) as io:
        io.add_primitive_shape("Disk", "Disk", [RADIUS])
        io.add_Newton_impact_friction_nsl("contact", mu=0.3, e=0)
        io.add_objects(
            "grains",
            [Contactor("Disk")],
            np.column_stack([x, np.zeros(n)]),
            velocities=velocities,
            mass=1.0,
        )

    options = sn.SolverOptions(sn.solver_ids.SICONOS_FRICTION_2D_NSGS)
    vnative_options = nonos.bridge.SpaceFilterOptions()
    vnative_options.neighborhood_radius = 3 * RADIUS
    vnative_options.min_radius = 0.5 * RADIUS

    run_options = MechanicsHdf5Runner_run_options()
    run_options["t0"] = 0
    run_options["T"] = 0.05
    run_options["h"] = 0.005
    run_options["theta"] = 0.5
    run_options["solver_options"] = options
    run_options["vnative_options"] = vnative_options
    run_options["verbose"] = False

    with MechanicsHdf5Runner(mode="r+", config=config, io_filename=filename) as io:
        io.run(run_options)
        work = np.asarray(io.contact_work_data())
        info = np.asarray(io.contact_info_data())

    assert work.shape[0] > 0
    # active interactions only, numbered as in cf_info: column 1 of
    # cf_work is column 1 of cf_info, which gives the body ids
    assert work.shape == (info.shape[0], 8)
    assert work[:, :2] == pytest.approx(info[:, :2])
    # the impact is seen: the fourth grain is in a contact at the end
    last = info[info[:, 0] == info[-1, 0]]
    assert n in last[:, 2:4]

    # theta = 1/2: the theta-weighted work is the mean velocity work
    assert work[:, 4] == pytest.approx(work[:, 2], abs=1e-14)
    assert work[:, 5] == pytest.approx(work[:, 3], abs=1e-14)
    # e = 0 and no contact separating before the step: no work of the
    # impact law
    assert np.all(work[:, 7] == 0)
    # the energy balance along x, where no force acts: the plastic impact
    # of one grain at 1 m/s on three leaves the four at 1/4 m/s, the
    # contact work is the kinetic energy change 4 * 1/2 * (1/4)^2 - 1/2
    # (the work of the velocity after the solve, not of the free one)
    assert work[:, 2].sum() + work[:, 3].sum() == pytest.approx(-0.375, abs=1e-6)
