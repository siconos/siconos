# cxd vendored copy

Vendored from the private cxd repository (`/home/maurice/src/cxd`, master),
commit `0f81b81` ([doc] README: no "fixed point vector"), 2026-10-19.
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
