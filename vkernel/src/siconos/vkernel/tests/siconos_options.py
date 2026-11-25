# Siconos is a program dedicated to modeling, simulation and control
# of non smooth dynamical systems.
#
# SPDX-FileCopyrightText: INRIA
# SPDX-License-Identifier: Apache-2.0
#

# Check: the simulation options of the siconos command line for the quasi
# inelastic shock law (--scheme), the rolling resistance (--mu-r) and the SVG
# pictures (--svg-*), on the vnative backend, through the runner.
#
# 1. The command line: the siconos script of the build tree parses a command
#    line into the options it passes to nonos/simoptions.py.
# 2. The options applied as the siconos command does (json, then
#    simoptions.apply) to a disk (R = 1, m = 1, I = 1/2) falling 0.5 onto a
#    ground segment, against the driver that chooses the same:
#
#   cli-mji        --scheme mji and SVG options, e = 0
#   osi-mji        osi = MoreauJeanQuasiInelasticOSI, e = 0     same bits
#   cli-mu-r       --mu-r 0.05, a Newton law, vx = 2
#   scene-mu-r     a rolling law mu_r = 0.05, vx = 2            same bits
#   cli-mu-r-0     --mu-r 0, a rolling law mu_r = 0.05, vx = 2
#   newton         a Newton law, vx = 2                         same bits
#   e-positive     --scheme mji, e = 0.9: refused
#   rolling-gosi   --scheme gosi --mu-r 0.05: refused
#
#    and the run options of the files: osi the class of the script, scheme
#    the scheme that ran; the SVG pictures in their directory, titled with
#    the scheme that ran, also by the siconos_svg converter.
#
# Each case runs in its own process: a second simulation in the same
# process crashes the vnative backend.
#
# usage: python3 siconos_options.py
# exit status 0 if every check passes

import importlib.machinery
import importlib.util
import json
import os
import subprocess
import sys
import tempfile

import h5py
import numpy as np

radius, m, inertia = 1.0, 1.0, 0.5
tables = ("dynamic", "velocities", "cf", "energy_work")


def run(case, filename, options):
    import siconos.numerics as sn
    from siconos.io.mechanics_run import (
        MechanicsHdf5Runner,
        MechanicsHdf5Runner_run_options,
        RunnerConfig,
    )
    from siconos.mechanics.collision.tools import Contactor

    import nonos  # noqa: F401 registers the vnative backend
    import nonos.bridge

    rolling = case in ("cli-mu-r", "scene-mu-r", "cli-mu-r-0", "newton")
    config = RunnerConfig("vnative")
    with MechanicsHdf5Runner(config=config, io_filename=filename) as io:
        io.add_primitive_shape("Disk", "Disk", [radius])
        io.add_primitive_shape("Ground", "Segment", (-10, 0, 10, 0))
        io.add_object(
            "disk",
            [Contactor("Disk")],
            translation=[0, radius + 0.5],
            orientation=[0],
            velocity=[2.0 if rolling else 0.0, 0, 0],
            mass=m,
            inertia=inertia,
        )
        io.add_object("ground", [Contactor("Ground")], translation=[0, 0])
        if case in ("scene-mu-r", "cli-mu-r-0"):
            io.add_Newton_impact_rolling_friction_nsl("contact", mu=0.5, mu_r=0.05, e=0)
        else:
            e = 0.9 if case == "e-positive" else 0.0
            io.add_Newton_impact_friction_nsl("contact", mu=0.5, e=e)

    # what the siconos command does with its simulation options
    if options:
        from nonos import simoptions

        simoptions.apply(json.loads(json.dumps(options)))

    so = sn.SolverOptions(sn.solver_ids.SICONOS_FRICTION_2D_NSGS)
    so.iparam[sn.params.SICONOS_IPARAM_MAX_ITER] = 100
    so.dparam[sn.params.SICONOS_DPARAM_TOL] = 1e-10

    run_options = MechanicsHdf5Runner_run_options()
    run_options["t0"] = 0
    run_options["T"] = 0.6
    run_options["h"] = 0.005
    run_options["theta"] = 0.50001
    run_options["solver_options"] = so
    run_options["output_frequency"] = 5
    run_options["output_energy_work"] = True
    if case == "osi-mji":
        run_options["osi"] = nonos.bridge.MoreauJeanQuasiInelasticOSI
    with MechanicsHdf5Runner(mode="r+", config=config, io_filename=filename) as io:
        io.run(run_options)


if len(sys.argv) == 4:
    case, filename, options = sys.argv[1:]
    run(case, filename, json.loads(options))
    sys.exit(0)


failed = False


def check(ok, message):
    global failed
    print(("ok: " if ok else "FAILED: ") + message)
    failed = failed or not ok


# 1. the command line, parsed by the siconos script of the build tree
script = os.path.join(os.environ["SICONOS_BUILD_DIR"], "scripts", "siconos")
loader = importlib.machinery.SourceFileLoader("siconos_script", script)
siconos_script = importlib.util.module_from_spec(
    importlib.util.spec_from_loader(loader.name, loader)
)
loader.exec_module(siconos_script)


def options_of(command_line):
    parser = siconos_script.parse_command_line()
    return siconos_script.collect_sim_options(parser.parse_args(command_line))


tmpdir = tempfile.mkdtemp(prefix="siconos_options_")
svg_dir = os.path.join(tmpdir, "pictures", "run")  # created by the run
cli_mji = options_of(
    [
        "--scheme",
        "mji",
        "--svg-output",
        "0.2,0.5",
        "--svg-dir",
        svg_dir,
        "--svg-coulomb",
        "--svg-animate",
        "--no-svg-cones",
        "--svg-material-side=-1=+1",
        "--svg-from",
        "0.1",
        "driver.py",
    ]
)
check(
    cli_mji
    == {
        "svg_output": "0.2,0.5",
        "scheme": "MoreauJeanQuasiInelastic",
        "svg_options": {
            "directory": svg_dir,
            "cones": False,
            "material_sides": [[-1, 1]],
            "coulomb": True,
            "animate": True,
            "window": [0.1, None],
        },
    },
    f"--scheme mji and the SVG options: {cli_mji}",
)
for alias, name in (
    ("quasi-inelastic", "MoreauJeanQuasiInelastic"),
    ("at-inelastic", "AnitescuTasoraInelastic"),
):
    scheme = options_of(["--scheme", alias, "driver.py"]).get("scheme")
    check(scheme == name, f"--scheme {alias}: {scheme}")
mu_r = options_of(["--mu-r", "0", "driver.py"])
check(mu_r == {"mu_r": 0.0}, f"--mu-r 0: {mu_r}")
check(options_of(["driver.py"]) == {}, "no option, no simulation option")


# 2. the options applied to a run
def run_case(case, options):
    filename = os.path.join(tmpdir, case + ".hdf5")
    r = subprocess.run(
        [sys.executable, __file__, case, filename, json.dumps(options)],
        capture_output=True,
        text=True,
    )
    return r, filename


def read(filename):
    with h5py.File(filename, "r") as f:
        data = f["data"]
        run_options = json.loads(data["siconos_mechanics_run_options"].attrs["options"])
        return {t: data[t][:] for t in tables if t in data}, run_options


def same(case_a, file_a, case_b, file_b):
    a, _ = read(file_a)
    b, _ = read(file_b)
    for t in sorted(set(a) | set(b)):
        ok = (
            t in a
            and t in b
            and a[t].shape == b[t].shape
            and np.array_equal(a[t], b[t])
        )
        check(ok, f"{case_a} and {case_b}: the same {t} table")


results = {}
for case, options in (
    ("cli-mji", cli_mji),
    ("osi-mji", {}),
    ("cli-mu-r", {"mu_r": 0.05}),
    ("scene-mu-r", {}),
    ("cli-mu-r-0", {"mu_r": 0.0}),
    ("newton", {}),
):
    r, filename = run_case(case, options)
    check(r.returncode == 0, f"{case}: exit status {r.returncode}")
    if r.returncode != 0:
        print(r.stdout[-2000:], r.stderr[-2000:])
    results[case] = filename

if not failed:
    same("cli-mji", results["cli-mji"], "osi-mji", results["osi-mji"])
    same("cli-mu-r", results["cli-mu-r"], "scene-mu-r", results["scene-mu-r"])
    same("cli-mu-r-0", results["cli-mu-r-0"], "newton", results["newton"])

    # the run options: the class of the script, the scheme that ran
    for case, osi in (
        ("cli-mji", "MoreauJeanOSI"),
        ("osi-mji", "MoreauJeanQuasiInelasticOSI"),
    ):
        _, run_options = read(results[case])
        check(
            run_options.get("osi") == osi
            and run_options.get("scheme") == "MoreauJeanQuasiInelastic",
            f"{case}: osi {run_options.get('osi')}, scheme {run_options.get('scheme')}",
        )
    # the rolling columns of the cf table with mu_r > 0 only
    columns = {
        case: read(results[case])[0]["cf"].shape[1]
        for case in ("cli-mu-r", "cli-mu-r-0")
    }
    check(columns == {"cli-mu-r": 32, "cli-mu-r-0": 26}, f"cf columns {columns}")
    # the disk is slowed down by the rolling resistance
    v_end = {
        case: read(results[case])[0]["velocities"][-1, 2]
        for case in ("cli-mu-r", "newton")
    }
    check(v_end["cli-mu-r"] < v_end["newton"], f"vx at the end {v_end}")

    # the SVG pictures, in their directory, titled with the scheme that ran
    names = sorted(os.listdir(svg_dir)) if os.path.isdir(svg_dir) else []
    expected = [
        "cli-mji_coulomb_000000.200000.svg",
        "cli-mji_coulomb_000000.500000.svg",
        "cli-mji_svg_000000.200000.svg",
        "cli-mji_svg_000000.500000.svg",
        "cli-mji_svg_animation.svg",
    ]
    check(names == expected, f"the pictures of the run: {names}")
    if names == expected:
        with open(os.path.join(svg_dir, "cli-mji_svg_000000.500000.svg")) as f:
            check(
                "MoreauJeanQuasiInelastic" in f.read(),
                "the title of a picture of the run",
            )
    # the converter reads the scheme that ran from the file
    from siconos.io.svg import convert

    converted = convert(
        results["cli-mji"], times=[0.5], directory=os.path.join(tmpdir, "conv")
    )
    with open(converted[0]) as f:
        check(
            "MoreauJeanQuasiInelastic" in f.read(), "the title of a converted picture"
        )

# the refusals: the checks of an osi of the script
for case, options, message in (
    ("e-positive", {"scheme": "mji"}, "replaces the restitution"),
    ("rolling-gosi", {"scheme": "gosi", "mu_r": 0.05}, "rolling resistance (mu_r > 0)"),
):
    r, _ = run_case(case, options)
    error = [line for line in r.stderr.splitlines() if line.startswith("siconos:")]
    check(
        r.returncode == 2 and any(message in line for line in error),
        f"{case} refused: exit status {r.returncode}, {error}",
    )

sys.exit(1 if failed else 0)
