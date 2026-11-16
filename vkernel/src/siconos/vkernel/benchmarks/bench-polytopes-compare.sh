#!/bin/sh
# settle N grains with vnative, dump positions, then time K steps on both
# vnative and bullet loading that static configuration. Prints a report at
# the end including memory, contact-detection time, and computeOneStep time.
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

# Extract timing and memory metrics from the compare output line.
# Format: ... contact_detection=X ms compute_one_step=X ms per_step=X ms memory=X MB
extract() {
  echo "$1" | sed -n "s/.*${2}=\([0-9.]*\)${3}.*/\1/p"
}

vs=$(extract "$vout" contact_detection "ms")
vc=$(extract "$vout" compute_one_step "ms")
vp=$(extract "$vout" per_step "ms")
vm=$(extract "$vout" memory "MB")

bs=$(extract "$bout" contact_detection "ms")
bc=$(extract "$bout" compute_one_step "ms")
bp=$(extract "$bout" per_step "ms")
bm=$(extract "$bout" memory "MB")

cat <<EOF

=================== POLYTOPES COMPARE REPORT ===================
N=${N} grains, K=${K} steps, static settled geometry

                       contact_detect   compute_step   per_step   memory
vnative :              ${vs} ms       ${vc} ms       ${vp} ms   ${vm} MB
bullet  :              ${bs} ms       ${bc} ms       ${bp} ms   ${bm} MB
EOF
if [ -n "$vp" ] && [ -n "$bp" ] && [ "$bp" != "0" ]; then
  awk -v v="$vp" -v b="$bp" 'BEGIN { printf "speedup (vnative/bullet): %.2fx\n", v/b }'
fi
echo "==============================================================="
