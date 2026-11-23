# Siconos 2026.rc3 Release Notes

**Release Date:** November 2026
**Previous Release:** 4.5.0 (release candidates 2026.rc1, 2026.rc2)
**Commits:** about 300 since 2026.rc2

---

## Overview

Siconos 2026.rc3 is the third release candidate of the 2026 series. It is mainly a
consolidation release: the examples of `siconos-tutorials` were ported to the refactored
kernel and run on macOS arm64 and on Ubuntu 26.04 x86_64, and this porting fixed many
numerical schemes, memory leaks, uninitialized values and platform-dependent tests. The
release also completes the renaming of the Lagrangian and Newton-Euler relations, makes the
FEM plasticity formulation well posed, adds a Levenberg-Marquardt fallback to the
semi-smooth Newton solvers, defines a versioned HDF5 format for the io component and brings
the vKernel to an installable component with new time-stepping schemes and polytope contacts.

The notes of the previous candidates ([RELEASE_2026.rc2.md](RELEASE_2026.rc2.md),
[RELEASE_2026.rc1.md](RELEASE_2026.rc1.md)) describe the changes since 4.5.0. The rc2 notes
were written after the rc2 tag: the API renamings listed there (relations, numerics waves 1
and 2) are part of 2026.rc3.

---

## What's New Since 2026.rc2

### 1. Kernel: Relations

- **Descriptive names for the relation matrices** of `LagrangianR` and `NewtonEulerR`:
  `H` is the contact matrix (`ẏ = H v`, input `p = Hᵀλ`, `contactMatrix()`); for
  Newton-Euler, `G = ∇_q h` is the constraint Jacobian (`constraintJacobian()`) and
  `H = G T`. Contact relation classes are named after their output
  (`Lagrangian2dNormalVelocityR`, `NewtonEulerVelocityR`, ...). See Breaking Changes.
- The former functions remain as deprecated aliases (a warning at the call); the former
  protected reference members (`jacobianhOver_q_view_`, `H_NE_view_`, ...) were removed.
- `NewtonEulerR` family: equations documented, contact Jacobians factored, bug fixed in
  `NewtonEulerVelocityAngularVelocityR` (with a new cppunit test).
- `Lagrangian2d3DR` (now `Lagrangian2dVelocityAngularVelocityR`): the rolling row lost in an
  earlier refactoring is restored.

### 2. Kernel: Time-Stepping and Event-Driven Schemes

- **Projected Moreau-Jean schemes** (`TimeSteppingCombinedProjection`,
  `TimeSteppingDirectProjection`): the index set iterations follow Alg. 2 of the CMAME 2013
  paper; new option `ProjectionUpdate::from_origin` for the exact projection
  `q = a + G(q)ᵀμ`. Documentation in `docs/devel/PROJECTION_SCHEMES.md`.
- **SchatzmanPaoliOSI**: start of the scheme, impact law and stiffness term fixed
  (`docs/devel/SCHATZMAN_PAOLI_SCHEME.md`).
- **D1MinusLinearOSI**: index sets computed at the start of the step, impact rule and
  tolerance options, support of `LagrangianSparseDS` in the `LinearOSNS` blocks.
- **EventDriven**: the acceleration of the dynamical systems is computed at the end time of
  the step (Lsodar, Hem5).
- `MoreauJeanOSI::updateInput` computes the input on `indexSet0`.
- `BlockVector::operator=` copies the values into the existing blocks (wrong results with
  `EulerMoreauOSI` and `FirstOrderType2R` / `FirstOrderNonLinearR`).
- `GlobalFrictionContact`: friction coefficient 0 for `MohrCoulombPlasticityNSL`.

### 3. Mechanics: FEM and Plasticity

- **Well-posed FEM plasticity formulation** (`SolidLinearTIDS`, `StressLinearTIR`,
  `MohrCoulombPlasticityNSL`): the iteration matrix is `W = [M, hθBᵀ; −hθB, S]` with `S` the
  compliance. The former symmetric form gave a negative Delassus operator for the stress
  interactions (no physical solution, solver failures). New accessors `plasticRate()`
  (`ε̇ₚ = −Cᵀλ/h`) and `plasticDeformation()`. Documentation in
  `docs/devel/FEM_PLASTICITY_SCHEME.md`.
- 3D `TH4` stiffness and mass matrices fixed; contact point of `NodeFem1d2DR` fixed.
- Heap overflow fixed when reading GMSH2 physical names.
- `SiconosBulletCollisionManager`: leak of the sphere-sphere collision function fixed.

### 4. Numerics

- **`GFC3D_NSN_AC`**: Alart-Curnier is now the default formulation (Jean-Moreau stalled on
  FCLIB problems with a nearly singular Jacobian), with a **Levenberg-Marquardt fallback** on
  stalled iterations (`SICONOS_FRICTION_3D_NSN_LM_FALLBACK`, `..._LM_COEF`). The helper
  (`nonsmooth_newton_levenberg_marquardt.{h,c}`) is shared with `GFC2D_NSN_AC` and the
  `FC3D_NSN_AC/FB/NM` solvers, where it is off by default. Documentation in
  `docs/devel/GFC3D_NSN_LM_FALLBACK.md`.
- **Sparse gemv** (`CSparseMatrix_aaxpby`, `SBM_gemv`, `NM_gemv`): `y` is overwritten when
  `beta == 0` (BLAS convention). A NaN in a reused work vector used to survive and made the
  FCLIB `gfc3d` collection fail randomly.
- **`MLCP_ENUM`**: second pass with an SVD (minimum-norm) solve of the singular cases, used
  only when the plain enumeration fails (redundant equality constraints).
- Rolling friction 2D: exact one-contact `ENUM` local solver.
- Cohesive NSGS: compliance, relaxation, Anderson acceleration, final status fixed, new test
  collection.
- `GFC_NSN`: option to compute ρ from the diagonal of an approximated `W`, exposed in Python.
- Parallel FC2D solvers: first-touch and no-copy variants, equitable graph coloring, clean
  failure of the CUDA solver when the coloring fails.
- Contact-status helpers moved to `FrictionContactProblem`; generic `nsgs_utils` helpers
  replace the problem-specific NSGS vector helpers; experimental batched NSGS passes.
- Public API cleanup (waves 1 and 2): cryptic abbreviations expanded, lifecycle verbs
  standardized, deprecated aliases in `numerics/src/utils/deprecated_numerics_api.h`.
- Memory leaks fixed in `CSparseMatrix`, the `NumericsMatrix` pybind11 interface and the NSGS
  local problems with `NM_SPARSE` storage.
- Hard-coded `verbose = 1` removed.

### 5. IO

- **Versioned HDF5 format**: `format_v1.json` (description of the current files) and
  `format_v2.json` (schema). The file stores its format version
  (`version_major`, `version_minor`), the format used (`/format`) and the script that created
  it (`data/script/source`). A file of another format version is refused.
  `mechanics_hdf5.py` uses `format_hdf5.py`, which checks that the data follow the format.
- **VTKHDF export** (ParaView) of bodies, shapes, static objects, contacts and FEM meshes.
- `vview`: non-disk bodies drawn with their shape in glyph mode.
- `mechanics_run`: the OSNSP is not solved when `preCompute` returns false; static objects of
  an output written at once; FEM body defaults and closed chains of FEM contact segments.

### 6. vKernel

The vKernel is now an installable component (`WITH_vkernel=ON`): namespace `siconos::vkernel`,
headers installed under `siconos/vkernel`, an exported `Siconos::vkernel` target, a source
tree that mirrors the install tree (`vkernel/src/siconos/vkernel`), and Apache license
headers. Main additions (details in `vkernel/ChangeLog`):

- **Time-stepping schemes chosen at run time**: Moreau-Jean (split and global/GOSI),
  Moreau-Jean with the Frémond impact law, Anitescu-Tasora (projected Gauss-Seidel and APGD),
  explicit DEM with penalty springs; `siconos --scheme ...` on the command line.
- **Polytope contacts**: narrowphase on the vendored `cxd` library, polytope/polytope,
  polytope/disk, polytope/segment and polytope/mesh-edge relations, contact manifolds by face
  clipping, per-grain shapes and convex polygon aggregates.
- **Storage**: registry only (the compile-time database and `NONOS_COMPILE_TIME_DATABASE`
  were removed), `storage::permute`, `storage::remove` with reference remapping,
  `storage::kernels` and `storage::schedule` (time step as kernels, with timings and
  dataflow).
- Several contact laws by collision group, polydisperse radii, `cuNSearch` fixes
  (reproducible sorts, no neighbor cap).
- Tests on cppunit and as ctest cases (NSCD analytic cases, comparisons MJ/AT/DEM),
  benchmarks (dam break of disks and rocks).

### 7. Build, Platforms and Tests

- **CMake component options**: components are selected with the standard
  `WITH_<component>` options (`WITH_numerics`, `WITH_kernel`, `WITH_mechanisms`,
  `WITH_vkernel`, ...). The `COMPONENTS` list of the configuration files still works.
- **macOS**: `BLA_VENDOR` defaults to Accelerate (the new LAPACK interface of macOS 13.3);
  the cached BLAS/LAPACK results are reset when the vendor changes.
- **Ubuntu 26.04** (GCC 15, CMake 4.2, Python 3.14, numpy 2.5, OpenBLAS, Boost 1.90): build
  and test fixes, Docker image in `ci_gitlab/dockerfiles/ubuntu26.04`. GCC 16 portability of
  the vKernel.
- **mechanisms**: builds again with OpenCASCADE 7.9 (deprecation warnings removed).
- Compilation warnings removed (kernel, numerics, externals) and ASan fixes; uninitialized
  Eigen objects fixed (results that differed on Linux only).
- **Non-regression tests of the time-stepping schemes** (`TimeSteppingSchemesTest`:
  bouncing beads, rocking block, cube, circuits...), validated against an independent Python
  implementation. `compareRefFile` checks the shapes of data and reference.
- `ZOHTest`, `ObserverTest`, `TwistingTest`: relative tolerances, comparisons that do not
  depend on the platform (they failed on macOS arm64).
- Guix: only the git-tracked files of the checkout are built; `lpsolve` from Guix.

---

## Breaking Changes Since 2026.rc2

1. **Lagrangian and Newton-Euler relations**: classes and functions renamed (table in
   [RELEASE_2026.rc2.md](RELEASE_2026.rc2.md), Breaking Changes, item 7). The former names
   compile with deprecation warnings; the former protected reference members were removed
   (compile error in derived classes that used them).
2. **Numerics function names**: waves 1 and 2 of the API cleanup (items 7 and 8 of the rc2
   notes); deprecated aliases in `deprecated_numerics_api.h`.
3. **`GFC3D_NSN_AC` default formulation** is Alart-Curnier (was Jean-Moreau): results and
   iteration counts of this solver change.
4. **Sparse gemv with `beta == 0`** overwrites `y` instead of scaling it.
5. **FEM plasticity**: the sign convention of the stress interactions changed with the
   iteration matrix; results of plasticity models change (they were wrong).
6. **HDF5 files**: a file of another format version is refused by `mechanics_hdf5`.
7. **vKernel**: namespace `siconos::vkernel`, include path `siconos/vkernel/...`,
   `NONOS_COMPILE_TIME_DATABASE` removed.
8. **CMake**: components with `WITH_<component>`; on macOS the default BLAS is Accelerate.

---

## Known Issues

- Python bindings of the `mechanisms` component (MBTB/CADMBTB) are not implemented: the
  `Mechanisms/*` tutorials cannot run.
- The OCC tutorials `occ_watch_escapement` and `occ_slider_crank` use C plugins loaded by
  name, which were removed from the kernel.
- `EventDriven`: the default tolerance (`DEFAULT_TOL_ED ≈ 2.2e-12`, absolute) is close to
  round-off; examples near Zeno accumulation set `setToleranceED(1e-8)`.
- Mass-weighted projection (`_useMassNormalization`) in the projected schemes is
  inconsistent.
- FEM plasticity: the non-global `fem::MoreauJeanOSI` path is untested, energy drift in
  `ballAgainstPlasticBernoulliBeam`.
- `lp_solve` 5.5 leaks on macOS (`realloc(p, 0)` in `LUSOL_free`).

---

## Testing

- macOS arm64 (Apple clang, Accelerate): all ctest tests pass (release build), ASan/UBSan
  build checked, leaks checked with `leaks --atExit`.
- Ubuntu 26.04 x86_64 (GCC 15, OpenBLAS): all ctest tests pass, except timeouts under
  emulation.
- `siconos-tutorials`: examples ported to the new API and run on both platforms, except the
  known issues above. The last fixes (`plasticBernoulliBeam`, RockFall
  `NSCD_model_adapted_energy_1702_FB`) were checked on macOS only.

---

## Contributors

Maurice Brémond, Vincent Acary, Michel Henry, Loris Megy, and the Siconos development team
at INRIA. See [credits.md](../../credits.md).

---

## Download & Installation

```bash
git clone https://gricad-gitlab.univ-grenoble-alpes.fr/nonsmooth/siconos.git
cd siconos
git checkout 2026.rc3

cmake -S . -B build -DCMAKE_BUILD_TYPE=Release
cmake --build build -j$(nproc)
cmake --install build
```

Use the `2026.rc3` tag of `siconos-tutorials` with this release.

---

## Documentation & Support

- **Documentation**: https://nonsmooth.gricad-pages.univ-grenoble-alpes.fr/siconos
- **Tutorials**: https://nonsmooth.gricad-pages.univ-grenoble-alpes.fr/siconos-tutorials
- **Issue Tracker**: https://gricad-gitlab.univ-grenoble-alpes.fr/nonsmooth/siconos/-/issues
- **Homepage**: http://siconos.org

---

**License**: Apache License 2.0
**Copyright**: 2024-2026 INRIA
