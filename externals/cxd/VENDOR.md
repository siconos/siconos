# cxd vendored copy

Vendored from the private cxd repository (`/home/maurice/src/cxd`, master),
commit `0f81b81` ([doc] README: no "fixed point vector"), 2026-10-19.

One local patch on top, not yet in the cxd remote (the user pushes it):
`src/cxd.hpp` + `src/cxd_batch.hpp` carry the 2D EPA degenerate-witness
bootstrap from the cxd branch `fix-epa-2d-degenerate-seed` (commit
`7066227`, worktree `/home/maurice/wkt/cxd/fix-epa-2d-degenerate-seed/cxd`):
on a zero-area GJK witness, EPA queries a support off the degenerate line
and restarts the face walk from that triangle, and the seed triangle's
orientation is canonicalized so the answer's normal does not depend on the
witness's vertex order. cxd's own suite is green with it (60/60, GCC and
Clang). Re-vendoring from master before that branch is pushed would drop
the fix and bring back the 2D EPA's NaN answers for flush overlaps.
cxd is not distributed yet; until then this directory is updated by hand
(copy the sources) and never re-formatted by the pre-commit hooks
(the `.pre-commit-config.yaml` exclusion pattern `^externals/` covers it).

Contents:

- `src/*.hpp`, `src/cxd.cpp` + `src/CMakeLists.txt` — host library
  (convex distance, batch GJK/EPA: `cxd_batch.hpp`, `cxd.hpp`, ...)
- `src/device/` — CUDA device backend (opt-in, not built here yet)
- `src/boost/hana` — the only boost dependency, vendored along
- `COPYING` — Apache 2.0 (all cxd sources carry the license header)

Excluded on purpose: `docs/`, `src/test/` (11M), demos, benchmarks,
build trees, shared-library artifacts.

Do not commit changes originated in this directory back to it; the cxd
upstream lives in `/home/maurice/src/cxd` and is never pushed to its
remote by agents (the user pushes himself).
