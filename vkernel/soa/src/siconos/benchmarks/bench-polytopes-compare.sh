#!/bin/sh
# settle N grains with vnative, dump positions, then time K steps on both
# vnative and bullet loading that static configuration. Prints a report at
# the end.
# Usage: bench-polytopes-compare.sh <N> [K]
[ -f bench-polytopes.py ] || { cd "$(dirname "$0")" || exit 1; }

N="${1?usage: bench-polytopes-compare.sh <N> [K]}"
K="${2:-20}"
PY="${PYTHON:-python3.11}"

# The narrowphase is serial. CompactNSearch runs its hash and query loops in
# OpenMP regions: at N=30 the 24-thread fork/join costs more than it saves
# (search 0.20 ms vs 0.06 ms with one thread). Bullet runs with one thread
# too, so both backends are compared single-threaded.
export OMP_NUM_THREADS=1

dump="settled-${N}.npz"

echo "### phase 1: settle vnative N=$N"
$PY bench-polytopes.py vnative "$N" --dump "$dump"

echo "### phase 2: compare $K steps on settled geometry"
vout=$($PY bench-polytopes.py vnative "$N" --load "$dump" --compare-steps "$K" 2>&1)
echo "$vout"
bout=$($PY bench-polytopes.py bullet "$N" --load "$dump" --compare-steps "$K" 2>&1)
echo "$bout"

# Extract per_step ms (value after "per_step=" up to the trailing "ms"): the
# mean of updateInteractions + computeOneStep over the timed loop.
vstep=$(echo "$vout" | sed -n 's/.*per_step=\([0-9.]*\)ms.*/\1/p')
bstep=$(echo "$bout" | sed -n 's/.*per_step=\([0-9.]*\)ms.*/\1/p')

cat <<EOF

=================== POLYTOPES COMPARE REPORT ===================
N=${N} grains, K=${K} steps, static settled geometry
vnative : ${vstep} ms/step
bullet  : ${bstep} ms/step
EOF
if [ -n "$vstep" ] && [ -n "$bstep" ] && [ "$bstep" != "0" ]; then
  awk -v v="$vstep" -v b="$bstep" 'BEGIN { printf "speedup (vnative/bullet): %.2fx\n", v/b }'
fi
echo "================================================================="
