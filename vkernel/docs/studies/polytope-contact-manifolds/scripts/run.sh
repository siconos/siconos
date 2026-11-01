#!/bin/sh
# The polytope contact manifold study: the scenes, the check against the
# simulator and the figures, with the vnative backend of a siconos build,
# in the guix container of guix.scm.
#
#   SICONOS_BUILD_DIR=<siconos build directory> ./run.sh [command...]
#
# Without a command: runs the seven scenes, check.py and figs.py. The
# outputs (hdf5, logs, png) go to $WORKDIR (default ./work).

set -e
here=$(cd "$(dirname "$0")" && pwd)
src=$(git -C "$here" rev-parse --show-toplevel)
build=$(cd "${SICONOS_BUILD_DIR:?set SICONOS_BUILD_DIR to a siconos build directory}" && pwd)
work=${WORKDIR:-$PWD/work}
mkdir -p "$work/pkg/nonos"
work=$(cd "$work" && pwd)

# the nonos package as the vnative runner imports it: the python sources
# and the extension module of the build
ln -sf "$src/vkernel/soa/src/siconos/py/nonos/__init__.py" \
  "$src/vkernel/soa/src/siconos/py/nonos/bridge.py" "$work/pkg/nonos/"
ln -sf "$build"/vkernel/soa/_nonos*.so "$work/pkg/nonos/"

if [ $# -eq 0 ]; then
  # shellcheck disable=SC2016 # expanded by the inner shell
  set -- sh -c 'for s in flat tilted overhang tipping junction stack corner; do
      echo "== $s"; python3 "$SCRIPTS/scenes.py" $s > $s.log 2>&1; done
    python3 "$SCRIPTS/check.py" flat tilted overhang tipping junction stack corner
    python3 "$SCRIPTS/figs.py"'
fi

# shellcheck disable=SC2016 # expanded inside the container
exec guix shell --container --share="$src" --share="$build" --share="$work" \
  -D -f "$src/guix.scm" python-h5py python-matplotlib bash coreutils -- \
  bash -c '
    b=$1; w=$2; export SCRIPTS=$3; shift 3
    libs=
    for d in control externals io kernel mechanics numerics vkernel/soa \
             vkernel/soa/extern/CompactNSearch; do
      libs=$libs$b/$d:
    done
    export LD_LIBRARY_PATH=$libs$GUIX_ENVIRONMENT/lib
    export PYTHONPATH=$b/python:$w/pkg:$SCRIPTS
    cd "$w"
    exec "$@"' run "$build" "$work" "$here" "$@"
