# SPDX-FileCopyrightText: 2026 INRIA
#
# SPDX-License-Identifier: Apache-2.0
#
# Siconos is a program dedicated to modeling, simulation and control
# of non smooth dynamical systems.
#
# The figures of vkernel/docs/notes/svg-display.md: the scenes of the
# time-integrators note (incline, sliding square, column, shear) and a
# granular one, run with the vnative backend, then drawn with
# siconos.io.svg.
#
# usage: python3 examples.py <directory> [case ...]
#   the cases: incline_stick incline_slide sliding column shear bounce
#   dambreak (all by default); the HDF5 files and the SVG files are written
#   in <directory>, the figures of the note are copied from there.
# Each simulation runs in its own process (a second vnative simulation in a
# process crashes). The python path is the one of the ctests (nonos and
# siconos.io).

import math
import os
import subprocess
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
DAMBREAK = os.path.join(
    HERE,
    "..",
    "..",
    "..",
    "..",
    "src",
    "siconos",
    "vkernel",
    "benchmarks",
    "bench-dambreak.py",
)
CASES = (
    "incline_stick",
    "incline_slide",
    "sliding",
    "column",
    "shear",
    "bounce",
    "dambreak",
)
SQUARE = np.array([[-0.5, -0.5], [0.5, -0.5], [0.5, 0.5], [-0.5, 0.5]])


def scene(case, filename, scheme):
    """Build and run a scene in this process."""
    import siconos.numerics as sn
    from siconos import integrators
    from siconos.mechanics.collision.tools import Contactor
    from siconos.io.mechanics_run import (
        MechanicsHdf5Runner,
        MechanicsHdf5Runner_run_options,
        RunnerConfig,
    )
    import nonos  # noqa: F401 registers the vnative backend
    import nonos.bridge

    config = RunnerConfig("vnative")
    T, h, every = 1.0, 0.005, 4
    with MechanicsHdf5Runner(config=config, io_filename=filename) as io:
        if case in ("incline_stick", "incline_slide"):
            alpha, mu = (20, 0.5) if case == "incline_stick" else (30, 0.3)
            a = math.radians(alpha)
            c, s = math.cos(a), math.sin(a)
            io.add_primitive_shape("Incline", "Segment", (-3 * c, 3 * s, 3 * c, -3 * s))
            io.add_object("incline", [Contactor("Incline")], translation=[0, 0])
            io.add_convex_shape("Square", SQUARE)
            # on the incline, centred at its middle
            io.add_object(
                "square",
                [Contactor("Square")],
                translation=[0.5 * s - 1.5 * c, 0.5 * c + 1.5 * s],
                orientation=[-a],
                velocity=[0, 0, 0],
                mass=1,
                inertia=1 / 6,
            )
            io.add_Newton_impact_friction_nsl("contact", mu=mu, e=0)
        elif case == "sliding":
            io.add_primitive_shape("Ground", "Segment", (-1, 0, 4, 0))
            io.add_object("ground", [Contactor("Ground")], translation=[0, 0])
            io.add_convex_shape("Square", SQUARE)
            io.add_object(
                "square",
                [Contactor("Square")],
                translation=[0, 0.5],
                orientation=[0],
                velocity=[2, 0, 0],
                mass=1,
                inertia=1 / 6,
            )
            io.add_Newton_impact_friction_nsl("contact", mu=0.3, e=0)
            T = 0.7
        elif case == "column":
            # five disks of radius 0.5 above each other, gaps of 0.01
            io.add_primitive_shape("Ground", "Segment", (-2, 0, 2, 0))
            io.add_object("ground", [Contactor("Ground")], translation=[0, 0])
            io.add_primitive_shape("Disk", "Disk", [0.5])
            for k in range(5):
                io.add_object(
                    f"disk{k}",
                    [Contactor("Disk")],
                    translation=[0, 0.51 + 1.01 * k],
                    orientation=[0],
                    velocity=[0, 0, 0],
                    mass=1,
                    inertia=0.125,
                )
            io.add_Newton_impact_friction_nsl("contact", mu=0.3, e=0)
            T, every = 0.3, 1
        elif case == "shear":
            # a layer of touching fixed disks of radius 0.5, a disk of the
            # same radius thrown over it from a hollow, its rotation nearly
            # locked: it climbs and leaves the crests (dilatancy)
            io.add_primitive_shape("Disk", "Disk", [0.5])
            for k in range(8):
                io.add_object(f"fixed{k}", [Contactor("Disk")], translation=[k, 0])
            io.add_object(
                "disk",
                [Contactor("Disk")],
                translation=[0.5, math.sqrt(3) * 0.5 + 0.001],
                orientation=[0],
                velocity=[5, 0, 0],
                mass=1,
                inertia=1e4,
            )
            io.add_Newton_impact_friction_nsl("contact", mu=0.1, e=0)
            T, every = 1.4, 2
        elif case == "bounce":
            io.add_primitive_shape("Floor", "Segment", (-1, 0, 16, 0))
            io.add_object("floor", [Contactor("Floor")], translation=[0, 0])
            io.add_primitive_shape("Disk", "Disk", [0.5])
            io.add_object(
                "disk",
                [Contactor("Disk")],
                translation=[0, 1.5],
                orientation=[0],
                velocity=[6, 0, 0],
                mass=1,
                inertia=0.125,
            )
            io.add_Newton_impact_friction_nsl("contact", mu=0.1, e=0.5)
            T, every = 2.6, 1

    options = sn.SolverOptions(sn.solver_ids.SICONOS_FRICTION_2D_NSGS)
    options.iparam[sn.params.SICONOS_IPARAM_MAX_ITER] = 10000
    options.dparam[sn.params.SICONOS_DPARAM_TOL] = 1e-10
    run_options = MechanicsHdf5Runner_run_options()
    run_options["t0"] = 0
    run_options["T"] = T
    run_options["h"] = h
    run_options["solver_options"] = options
    run_options["output_frequency"] = every
    run_options["verbose"] = False
    run_options["osi"] = {
        "mj": integrators.MoreauJeanOSI,
        "at": nonos.bridge.AnitescuTasoraOSI,
    }[scheme]
    with MechanicsHdf5Runner(mode="r+", config=config, io_filename=filename) as io:
        io.run(run_options)


def simulate(directory, case, scheme="mj"):
    filename = os.path.join(directory, f"{case}_{scheme}.hdf5")
    subprocess.run(
        [sys.executable, __file__, "--run", case, filename, scheme],
        check=True,
        cwd=directory,
    )
    return filename


def draw(directory, case):
    """The figures of a case: returns their file names."""
    from siconos.io import svg

    def convert(filename, name, **options):
        names = svg.convert(filename, directory=directory, **options)
        out = []
        for old in names:
            kind = os.path.basename(old)[
                len(os.path.splitext(os.path.basename(filename))[0]) :
            ]
            new = os.path.join(directory, f"{name}{kind}")
            os.replace(old, new)
            out.append(new)
        return out

    if case == "incline_stick":
        f = simulate(directory, case)
        return convert(f, case, cones=True)
    if case == "incline_slide":
        f = simulate(directory, case)
        return convert(f, case, cones=True, animate=True, strobe=True, stride=5)
    if case == "sliding":
        f = simulate(directory, case)
        return convert(f, case, strobe=True, stride=4, times=[0.3], cones=True)
    if case == "column":
        f = simulate(directory, case)
        return convert(f, case, animate=True, uniform=True, discrete=True, times=[0.3])
    if case == "shear":
        f = simulate(directory, case)
        return convert(f, case, strobe=True, stride=4, animate=True, colour="none")
    if case == "bounce":
        f = simulate(directory, case)
        # the first impact, at t = 0.455
        return convert(f, case, animate=True, times=[0.455], coulomb=True)
    if case == "dambreak":
        # the scene of bench-dambreak.py at size 40, its own times: the
        # grains settle behind the gate, the gate is removed, they flow
        size = 40
        f = os.path.join(directory, f"bench-dambreak-{size}.hdf5")
        subprocess.run(
            [sys.executable, os.path.normpath(DAMBREAK), str(size)],
            check=True,
            cwd=directory,
        )
        import h5py

        with h5py.File(f, "r") as io:
            T = float(io["data/step_index/dynamic"][-1, 0])
        names = convert(
            f,
            "dambreak",
            times=[T / 3, 2 * T / 3, T],
            coulomb=True,
            chains=True,
            arrows=False,
        )
        # the front on the slope, cropped (said in the title)
        names += convert(
            f,
            "dambreak_front",
            times=[2 * T / 3],
            view=(420, 720, -30, 110),
            arrows=False,
            chains=True,
        )
        names += convert(f, "dambreak", animate=True, stride=12, colour="none")
        return names
    raise ValueError(case)


if __name__ == "__main__":
    if sys.argv[1] == "--run":
        scene(sys.argv[2], sys.argv[3], sys.argv[4])
        sys.exit(0)
    directory = os.path.abspath(sys.argv[1])
    os.makedirs(directory, exist_ok=True)
    for case in sys.argv[2:] or CASES:
        for name in draw(directory, case):
            print(f"{name}: {os.path.getsize(name)} bytes")
