#!/bin/sh
# run.py in the guix container of guix.scm, with matplotlib for the figures:
#
#   SICONOS_BUILD_DIR=<siconos build directory> ./run.sh [scene...]
#
# The outputs (hdf5, logs, png) go to $WORKDIR (default ./work). In the
# environment of the build, the cmake target
# vkernel_study_polytope_contact_manifolds does the same.

set -e
here=$(cd "$(dirname "$0")" && pwd)
src=$(git -C "$here" rev-parse --show-toplevel)
build=$(cd "${SICONOS_BUILD_DIR:?set SICONOS_BUILD_DIR to a siconos build directory}" && pwd)
work=${WORKDIR:-$PWD/work}
mkdir -p "$work"
work=$(cd "$work" && pwd)

exec guix shell --container --share="$src" --share="$build" --share="$work" \
  -D -f "$src/guix.scm" python-h5py python-matplotlib -- \
  python3 "$here/run.py" --build-dir "$build" --out "$work" "$@"
