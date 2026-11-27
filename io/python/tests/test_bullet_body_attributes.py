# SPDX-FileCopyrightText: INRIA
#
# SPDX-License-Identifier: Apache-2.0
#
# Siconos is a program dedicated to modeling, simulation and control
# of non smooth dynamical systems.
"""The attributes of a bullet body given to add_object are applied by the runner:
fext replaces the constant external force (the weight), as setConstantFext;
allow_self_collide: the bodies jointed to a body that does not allow it do not
collide with it."""

import numpy as np
import pytest

has_bullet = False
try:
    import siconos.io.mechanics_run as smrun
    import siconos.mechanics.collision.bullet  # noqa: F401
    import siconos.mechanics.collision.tools as smct

    has_bullet = True
except ImportError:
    pass

pytestmark = pytest.mark.skipif(not has_bullet, reason="Only if bullet is available")


def run(filename, build, T=0.2, h=5e-3):
    """The scene written by build, run; its dynamic, velocities and cf tables."""
    config = smrun.RunnerConfig(backend="bullet")
    with smrun.MechanicsHdf5Runner(io_filename=filename, config=config) as io:
        build(io)
    options = smrun.MechanicsHdf5Runner_run_options()
    options["T"], options["h"], options["verbose"] = T, h, False
    with smrun.MechanicsHdf5Runner(
        io_filename=filename, mode="r+", config=config
    ) as io:
        io.run(options)
        cf = io.contact_forces_data()
        return (
            io.dynamic_data()[()],
            io.velocities_data()[()],
            None if cf is None else cf[()],
        )


def ball(fext):
    def build(io):
        io.add_primitive_shape("Ball", "Sphere", (0.1,))
        io.add_object(
            "ball",
            [smct.Contactor("Ball")],
            translation=[0, 0, 1],
            mass=1,
            inertia=[1, 1, 1],
            fext=fext,
        )

    return build


def test_fext_replaces_the_weight_3d(tmp_path):
    # fext = 0: no force at all, the ball does not move
    dynamic, _, _ = run(str(tmp_path / "zero.hdf5"), ball([0, 0, 0]))
    assert np.abs(dynamic[:, 2:5] - [0, 0, 1]).max() < 1e-9
    # fext = (1, 0, 0), mass 1: x(t) = t^2 / 2, z constant
    dynamic, _, _ = run(str(tmp_path / "push.hdf5"), ball([1, 0, 0]))
    t = dynamic[-1, 0]
    assert dynamic[-1, 2] == pytest.approx(t * t / 2, abs=2e-3)
    assert np.abs(dynamic[:, 4] - 1).max() < 1e-9


def test_fext_2d(tmp_path):
    # fext = (fx, fy, torque) = (0, 0, 0.5), inertia 0.5: y constant, w(t) = t
    def build(io):
        io.add_primitive_shape("Disk", "Disk", (0.1,))
        io.add_object(
            "disk",
            [smct.Contactor("Disk")],
            translation=[0, 1],
            orientation=[0],
            velocity=[0, 0, 0],
            mass=1,
            inertia=0.5,
            fext=[0, 0, 0.5],
        )

    dynamic, velocities, _ = run(str(tmp_path / "disk.hdf5"), build)
    assert np.abs(dynamic[:, 3] - 1).max() < 1e-9
    assert velocities[-1, -1] == pytest.approx(dynamic[-1, 0], abs=1e-2)


def pair(allow_self_collide):
    """Two balls in contact, jointed by a pivot that allows self collision."""

    def build(io):
        io.add_primitive_shape("Ball", "Sphere", (0.1,))
        io.add_Newton_impact_friction_nsl("contact", mu=0.3, e=0.0)
        kw = (
            {}
            if allow_self_collide is None
            else {"allow_self_collide": allow_self_collide}
        )
        io.add_object(
            "a", [smct.Contactor("Ball")], translation=[0, 0, 0], mass=1, **kw
        )
        io.add_object(
            "b", [smct.Contactor("Ball")], translation=[0.15, 0, 0], mass=1, **kw
        )
        io.add_joint(
            "pivot",
            "a",
            "b",
            [[0.075, 0, 0]],
            [[0, 1, 0]],
            "PivotJointR",
            absolute=True,
            allow_self_collide=True,
        )

    return build


@pytest.mark.parametrize("allow_self_collide, collide", [(None, True), (False, False)])
def test_allow_self_collide(tmp_path, allow_self_collide, collide):
    # default (True, as in C++): the jointed balls collide; False: they do not
    _, _, cf = run(str(tmp_path / "pair.hdf5"), pair(allow_self_collide), T=0.05)
    contacts = 0 if cf is None else cf.shape[0]
    assert (contacts > 0) == collide
