# Siconos 2026-11 Release Notes

**Release Date:** November 2026
**Tag:** `2026-11` (version number `2026.11` in the `VERSION` file)
**Previous Release:** 4.5.0 (March 2024)
**Release candidates:** 2026.rc1, 2026.rc2, 2026.rc3
**Commits:** about 2000 since 4.5.0 (about 110 since 2026.rc3)

---

## Overview

Siconos 2026-11 is the official release of the 2026 series. It is a major release: the C++
kernel was modernized (C++20, `siconos::` namespaces, Eigen algebra), the Python bindings were
rewritten with pybind11, the numerics solvers were reorganized around a registration system,
and a new experimental component, the vKernel, simulates large-scale systems with a
structure-of-arrays storage. The release also brings cohesive zone models, plasticity solvers,
rolling friction, parallel and GPU solvers, a versioned HDF5 format, and a consolidation of the
time-stepping schemes checked by non-regression tests on macOS arm64 and Linux x86_64.

Siconos now uses calendar versioning: the version number (`YEAR.MONTH[.PATCH]`) is written
in the `VERSION` file at the root of the sources, read by CMake (CMake package
`siconos-2026.11`) and copied to `codemeta.json`.

The details are in the notes of the release candidates:

- [docs/devel/RELEASE_2026.rc1.md](docs/devel/RELEASE_2026.rc1.md): vKernel, solver
  registration, pybind11, C++20 modernization;
- [docs/devel/RELEASE_2026.rc2.md](docs/devel/RELEASE_2026.rc2.md): cohesive zone models,
  parallel and GPU solvers, GFC2D solvers, API renaming tables;
- [docs/devel/RELEASE_2026.rc3.md](docs/devel/RELEASE_2026.rc3.md): relations renaming,
  fixes of the time-stepping schemes, FEM plasticity, Levenberg-Marquardt fallback, versioned
  HDF5 format, installable vKernel.

---

## Highlights Since 4.5.0

### Kernel (C++)
- C++20, classes in the `siconos::` namespaces, `std::shared_ptr` instead of `SP::`.
- Eigen-based algebra (`SiconosMatrix`, `SiconosVector`, `BlockVector`, sparse matrices).
- User functions as `std::function` (lambdas, Python callables) instead of plugins loaded by
  name.
- Lagrangian and Newton-Euler relations with descriptive names (`contactMatrix`,
  `constraintJacobian`, `Lagrangian2dVelocityR`, `NewtonEulerVelocityR`, ...).
- Time-stepping schemes fixed and checked by non-regression tests: projected Moreau-Jean
  (exact projection option), Schatzman-Paoli, D1MinusLinear, event-driven, Lagrangian sparse
  relations, boundary conditions.

### Numerics (C)
- Solver registration (`REGISTER_SOLVER`), unified drivers, public API cleanup.
- Rolling friction (2D/3D, global), GFC2D solvers, cohesive friction contact, plasticity
  (Von Mises, Drucker-Prager) and mixed friction/plasticity problems.
- Interior point methods, semi-smooth Newton with a Levenberg-Marquardt fallback
  (`GFC3D_NSN_AC`, Alart-Curnier by default).
- Parallel (OpenMP, graph coloring, PETSc) and GPU (CUDA) solvers, optional.

### Mechanics and IO
- Cohesive zone models, FEM (native T3/TH4 elements, contact, plasticity).
- FEM plasticity with Von Mises (plane stress, plane strain) and Drucker-Prager
  (Mohr-Coulomb in the plane) laws, solved by a global semi-smooth Newton method and
  validated against pytripopfem (strip with a hole, footing).
- Versioned HDF5 format (`format_v1.json`, `format_v2.json`), VTKHDF export for ParaView,
  `siconos_filter`, `siconos_svg`, `vview` improvements.

### Python
- pybind11 bindings for numerics, kernel, mechanics, control and io (SWIG removed),
  Eigen/NumPy mapping with explicit `copy_t`/`alias_t` semantics.

### vKernel (experimental)
- Structure-of-arrays storage, run-time schemes (Moreau-Jean, Frémond, Anitescu-Tasora,
  explicit DEM, Jean quasi-inelastic law), disks, polytopes and FEM bodies, rolling friction,
  installable component (`WITH_vkernel=ON`).

### Build
- `WITH_<component>` CMake options, scikit-build for plugins, Accelerate by default on macOS,
  Ubuntu 26.04 support (GCC 15, Python 3.14, numpy 2.5), SPDX license headers and REUSE
  compliance, pre-commit hooks (`prek`).

---

## What's New Since 2026.rc3

### Kernel
- **Plasticity laws**: `VonMisesPlasticityNSL(sigma_y, hypothesis)` (J2, plane stress or
  plane strain) and `DruckerPragerPlasticityNSL(eta, theta, cohesion)` with
  `mohrCoulombPlaneStrain(c, phi, psi)`, in Mandel notation (documented in the headers).
- `OSNSMatrix::fillH` no longer zeroes the stress columns of a `StressLinearTIR` at the
  indices of the prescribed velocities (stress interactions are not concerned by velocity
  boundary conditions).
- Sparse iteration matrix property of the dynamical systems
  (`DynamicalSystemProperties::iterationMatrixSparse`), used by `OSNSMatrix::fillW` for the
  global one-step problems.
- `SiconosGraph`: vertices and edges renumbered only after a structural change (the
  renumbering is O(V + E) and was done at each step and each Newton iteration; FEM elements
  sharing one dynamical system give O(n²) edges).

### Mechanics: FEM plasticity
Documentation in `docs/devel/FEM_PLASTICITY_SCHEME.md`.
- `Material(rho, E, nu, AnalysisType2D, thickness)`: plane stress or plane strain, thickness;
  the compliance `S = D⁻¹` follows the analysis type (it was always plane strain);
  lumped mass (`MassType::lumped`).
- `fem::nonsmooth_formulations::Plasticity`: one-step problem of the FEM plasticity, condensed
  on the constrained stresses (built once at the first step), solved by the numerics
  plasticity driver with `PLASTICITY_2D_NSN` by default; Von Mises or Drucker-Prager laws;
  `setProblemDump` writes the numerics problems of chosen steps.
- `SolidLinearTIDS`: stress and strain unknowns in Mandel notation
  (`sxx, syy, √2 sxy`); `elementStressSelection` gives the physical stress of an element.
- `fem::MoreauJeanGOSI`: non-zero prescribed velocities fixed (they grew at every step with
  the wrong sign), sparse iteration matrix, no LU of the iteration matrix (unused).
  Memory and time on the strip with 1446 T3: 1.14 → 0.92 GB, 167 s → 104 s.
- Validated against pytripopfem: strip with a circular hole (Von Mises, plane stress) and
  half footing (Drucker-Prager, associated and non-associated flow) to round-off.

### Numerics
- **`PLASTICITY_2D_NSN`**: global semi-smooth Newton on the natural map of the whole
  plasticity problem (generalized Jacobian with the sparsity of M, sparse LU), full Newton
  steps with a watchdog by default, merit-function globalization as an option. NSGS needs
  O(cond) sweeps in quasi-static regimes.
- Plasticity `.dat` files keep the stress coordinates, the Von Mises hypothesis and the
  cohesion; matrices written with 17 digits; FEM strip and footing problems as test data.
- **Plasticity** (`numerics/src/Plasticity`, documentation in `docs/devel/PLASTICITY.md`):
  2D plasticity in Cartesian stress coordinates, Von Mises (J2, including the plane stress
  model) and Drucker-Prager criteria with cohesion, local solvers by return mapping and
  semi-smooth Newton, relative local errors by default, Python driver, `compute_error` and
  solve tests.
- **MixedProblems** (`numerics/src/MixedProblems`): problems mixing friction contacts and
  plasticity elements, element-wise NSGS on the global matrix.
- Friction solvers keep the options given by the user; parallel FC2D solvers fixed.
- Rolling friction 2D `ENUM`: fast path for a block without rolling row, rejection of a
  closed-form case that is not finite.
- The short-name macros of the friction solver IDs (`FC3D_NSGS`, `OC_NSN`, ...) were removed:
  use the `SICONOS_FRICTION_*` names.

### IO
- `siconos_filter`: copy of v1 files (scalar datasets, attributes, empty tables, numpy 2),
  time filter of the step-indexed tables, `--single` converts the floats only, `--exclude`
  also removes the joints, permanent interactions, boundary conditions and `cf_info` rows of
  an excluded object.
- `siconos_svg`: SVG pictures and animations of 2D scenes.
- `MechanicsIO` (`contactPoints`, `contactInfo`, `contactContactWork`,
  `contactInternalVariables`): the width of the contact rows comes from the contacts. With
  joints and contacts in the same index set, the rows were truncated to 0 columns and the
  time was written in every column of `cf`. `contactContactWork` accepts `(gap_law, h)` with
  every backend.

### vKernel
- Rolling friction for disks (Moreau-Jean, Moreau-Jean-Frémond, Anitescu-Tasora) and
  polytopes, checked against the kernel.
- Jean quasi-inelastic shock law with its Moreau-Jean local scheme; signed gap
  stabilisation for Anitescu-Tasora inelastic.
- FEM Dirichlet conditions with the free-dof operator, cold restart of the FEM bodies,
  energy/work outputs.

### Build and tests
- Calendar versioning (`VERSION` file).
- SPDX headers instead of the Apache-2.0 boilerplate, REUSE lint at push; `prek` runs the
  formatting and lint hooks at commit.
- Optional validation of the time-stepping schemes against hots-python
  (`-DWITH_HOTS_VALIDATION=ON`); configuration file with all `WITH_*` options.

---

## Breaking Changes Since 4.5.0

The migration tables are in the rc notes; in short:

1. C++: namespaces (`siconos::...`), `std::shared_ptr` instead of `SP::`, Eigen matrices
   (not initialized: call `setZero()`), user functions as `std::function`.
2. Relations renamed (rc2 notes, Breaking Changes, item 7); the former function names
   compile with deprecation warnings.
3. Numerics: function and solver renamings (rc2 notes, items 7 and 8; aliases in
   `deprecated_numerics_api.h`), `Plasticity2DProblem` instead of `MohrCoulomb2DProblem`,
   short-name macros of the friction solver IDs removed.
4. Python: pybind11 bindings, new module layout (`siconos.modeling`, `siconos.simulation`,
   `siconos.integrators`, `siconos.nonsmooth_formulations`, ...), `copy_t`/`alias_t`
   arguments, matrices in column-major order.
5. `GFC3D_NSN_AC` uses Alart-Curnier by default; sparse gemv overwrites `y` when
   `beta == 0`; FEM plasticity sign convention changed.
6. FEM: `SolidLinearTIDS` stress and strain unknowns in Mandel notation (`√2 sxy`); results
   with non-zero prescribed velocities in `fem::MoreauJeanGOSI` change (they were wrong).
7. HDF5: a file of another format version is refused by `mechanics_hdf5`.
8. CMake: `WITH_<component>` options, version `2026.11` (CMake package `siconos-2026.11`).

---

## Known Issues

- **siconos-tutorials, Python control examples** (`examples/control/*/*.py`): not updated
  to the pybind11 bindings (former `siconos.kernel`/`siconos.control` modules, former
  constructors, missing `VI` and `MixedComplementarityProblem2` numerics bindings), and no
  control example is run by the tutorials tests. The C++ control examples run, but some
  references are old (`TwoLinkManipulator_{Flexible, MultiConstraints,
  MultiConstraintsFlexible}` differ from their 2018 references).
- **siconos-tutorials, Jupyter notebooks**: not updated (`siconos-notebooks/st01`–`st03`
  call constructors without `copy_t`/`alias_t`, `advanced_01_impacting_bar` and
  `DiodeBridge.ipynb` import `siconos.kernel`, `advanced_02_guitar_string` needs
  `numpywrappers`). The Python scripts of the same examples are up to date.
- `mechanisms`: no Python bindings (MBTB/CADMBTB); the OCC examples using C plugins loaded
  by name (`occ_watch_escapement`, `occ_slider_crank`) cannot run.
- `EventDriven`: the default tolerance (`DEFAULT_TOL_ED ≈ 2.2e-12`, absolute) is close to
  round-off; examples near Zeno accumulation set `setToleranceED(1e-8)`.
- Mass-weighted projection (`_useMassNormalization`) in the projected schemes is
  inconsistent.
- FEM plasticity: on fine meshes the kernel dominates the cost (interaction graph with
  O(n²) edges, dense matrix `C` of `StressLinearTIR`); the relaxation option of NSGS is
  wrong for plasticity (applied after the projection); the non-global
  `fem::MoreauJeanOSI` path is untested.
- `lp_solve` 5.5 leaks on macOS (`realloc(p, 0)` in `LUSOL_free`).

---

## Testing

- macOS arm64 (Apple clang, Accelerate): all ctest tests pass (release build); ASan/UBSan
  build and leak checks done for 2026.rc3.
- Ubuntu 26.04 x86_64 (GCC 15, OpenBLAS): all ctest tests passed for 2026.rc3, except
  timeouts under emulation.
- siconos-tutorials 2026-11: the 332 tests (C++ and Python examples) pass on macOS arm64;
  FEM plasticity examples `Plasticity/StripWithHole` and `Plasticity/Footing` compared with
  pytripopfem.

---

## Contributors

Maurice Brémond, Vincent Acary, Franck Pérignon, Hoang-Minh Nguyen, Samuel Heidmann, Olivier
Goury, Loris Megy, Michel Henry, Paul Armand, Joël Chebassier, and the Siconos development
team at INRIA. See [AUTHORS](AUTHORS) and [credits.md](credits.md).

---

## Download & Installation

```bash
git clone https://gricad-gitlab.univ-grenoble-alpes.fr/nonsmooth/siconos.git
cd siconos
git checkout 2026-11

cmake -S . -B build -DCMAKE_BUILD_TYPE=Release
cmake --build build -j$(nproc)
cmake --install build
```

Use the `2026-11` tag of siconos-tutorials with this release.

---

## Documentation & Support

- **Documentation**: https://nonsmooth.gricad-pages.univ-grenoble-alpes.fr/siconos
- **Tutorials**: https://nonsmooth.gricad-pages.univ-grenoble-alpes.fr/siconos-tutorials
- **Issue Tracker**: https://gricad-gitlab.univ-grenoble-alpes.fr/nonsmooth/siconos/-/issues
- **Homepage**: http://siconos.org

---

**License**: Apache License 2.0
**Copyright**: INRIA
