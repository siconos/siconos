"""The polytope contact manifold study: the scenes (vnative, one process
each), the check against the simulator and the figures, written to the
output directory.

  python3 run.py --build-dir <siconos build> [--out <dir>] [scene...]

The cmake target vkernel_study_polytope_contact_manifolds runs it with the
paths of the build. The figures need matplotlib (in the development
environment of guix.scm): without it, they are skipped."""

import argparse
import glob
import importlib.util
import os
import subprocess
import sys

SCENES = ["flat", "tilted", "overhang", "tipping", "junction", "stack", "corner"]
LIBRARY_DIRS = [
    "control",
    "externals",
    "io",
    "kernel",
    "mechanics",
    "numerics",
    "vkernel/soa",
    "vkernel/soa/extern/CompactNSearch",
    "vkernel/soa/extern/cuNSearch",
]

here = os.path.dirname(os.path.abspath(__file__))
# the source tree: <src>/vkernel/docs/studies/<study>/scripts
src = os.path.abspath(os.path.join(here, *[".."] * 5))

parser = argparse.ArgumentParser(
    description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter
)
parser.add_argument("--build-dir", required=True, help="the siconos build directory")
parser.add_argument(
    "--out", default="work", help="the output directory (default: work)"
)
parser.add_argument(
    "--nonos-module", help="the _nonos extension (default: found in the build)"
)
parser.add_argument(
    "--library-path",
    default="",
    help="the directories of the siconos libraries, ':'-separated "
    "(default: the component directories of the build)",
)
parser.add_argument("scenes", nargs="*", default=SCENES)
args = parser.parse_args()

build = os.path.abspath(args.build_dir)
out = os.path.abspath(args.out)
nonos = (
    args.nonos_module
    or (glob.glob(os.path.join(build, "vkernel/soa/_nonos*.so")) or [None])[0]
)
if nonos is None:
    sys.exit(f"no _nonos module in {build}/vkernel/soa: build the target _nonos")

# the nonos package as the vnative runner imports it: the python sources
# and the extension module of the build
pkg = os.path.join(out, "pkg", "nonos")
os.makedirs(pkg, exist_ok=True)
for f in (
    os.path.join(src, "vkernel/soa/src/siconos/py/nonos/__init__.py"),
    os.path.join(src, "vkernel/soa/src/siconos/py/nonos/bridge.py"),
    os.path.abspath(nonos),
):
    link = os.path.join(pkg, os.path.basename(f))
    if os.path.lexists(link):
        os.remove(link)
    os.symlink(f, link)

libs = [d for d in args.library_path.split(":") if d] or [
    os.path.join(build, d) for d in LIBRARY_DIRS
]
# in a guix environment other than the one of the build, the RUNPATH of the
# libraries does not lead to the dependencies: the environment's lib does
if os.environ.get("GUIX_ENVIRONMENT"):
    libs.append(os.path.join(os.environ["GUIX_ENVIRONMENT"], "lib"))
env = dict(os.environ)
env["LD_LIBRARY_PATH"] = ":".join(libs + [p for p in [env.get("LD_LIBRARY_PATH")] if p])
env["PYTHONPATH"] = ":".join(
    [os.path.join(build, "python"), os.path.join(out, "pkg"), here]
    + [p for p in [env.get("PYTHONPATH")] if p]
)


def run(*cmd, log=None):
    print(" ".join(cmd), flush=True)
    if log:
        with open(log, "w") as f:
            r = subprocess.run(
                cmd, cwd=out, env=env, stdout=f, stderr=subprocess.STDOUT
            )
    else:
        r = subprocess.run(cmd, cwd=out, env=env)
    if r.returncode != 0:
        sys.exit(
            f"failed ({r.returncode}): {' '.join(cmd)}"
            + (f", see {log}" if log else "")
        )


# one process per scene: a second vnative simulation in the same process
# is not supported
for s in args.scenes:
    run(
        sys.executable,
        os.path.join(here, "scenes.py"),
        s,
        log=os.path.join(out, s + ".log"),
    )
run(sys.executable, os.path.join(here, "check.py"), *args.scenes)
if set(args.scenes) != set(SCENES):
    print("figures: they need all the scenes")
elif importlib.util.find_spec("matplotlib") is None:
    print("figures: skipped, no matplotlib")
else:
    run(sys.executable, os.path.join(here, "figs.py"))
    print(f"figures in {out}")
