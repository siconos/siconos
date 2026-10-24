# Siconos 2026.rc2 Release Notes

**Release Date:** October 2026
**Previous Release:** 4.5.0 (2026.rc1 was the first release candidate)
**Commits:** 1600+ changes since 4.5.0 (407 since 2026.rc1)

---

## Overview

Siconos 2026.rc2 is the second release candidate for the 2026 series. It builds on the major architectural work introduced in 2026.rc1 and adds significant new capabilities in numerics, vKernel, mechanics, and build quality. The release further modernizes the C++ and Python APIs, introduces cohesive zone modeling, experimental GPU and parallel CPU solvers, and hardens the solver and matrix infrastructure.

---

## Major Highlights Since 4.5.0

### 1. vKernel: New Architecture for Large-Scale Systems

Introduction of **vKernel** (virtual kernel), an experimental parallel architecture to the classic kernel, designed for high-performance simulation of large-scale nonsmooth systems:

- **SoA (Structure of Arrays)** data layout for cache-friendly, vectorized memory access
- **Runtime-defined systems** with flexible topology management (modern C++20)
- **Eigen-based algebra** backend with type-safe matrix/vector concepts
- **Registry-backed storage** with dynamic add/remove of objects and properties
- **FEM and disk contact support**: disk/fem interaction assembly, contact forces output, and `vview` visualization

**Key directories:**
- `vkernel/soa/src/siconos/algebra/`
- `vkernel/soa/src/siconos/storage/`
- `vkernel/soa/src/siconos/topology/`
- `vkernel/soa/src/siconos/simul/`

### 2. Numerics: Solver Registration & New Solvers

A complete refactoring of the solver infrastructure:

- **Solver registration system** (`solver_registry.h`) with `REGISTER_SOLVER` macros for runtime solver lookup
- **Unified driver architecture** across all solver families
- **Rolling friction contact solvers** (3D/2D)
- **Global rolling friction solvers** (IPM, NSGS-WR)
- **Plasticity 2D** module (renamed from Mohr-Coulomb)
- **GFC2D solvers**: new 2D global friction-contact solvers including NSN-AC, FPP, and EG variants
- **User-provided rho vectors** for `GFC3D_NSN_AC` and `fc3d_compute_rho` helper with Python bindings
- **Improved NSN-AC warm start, regularization, and diagnostics**
- **Standard verbose mode** and better sparsity reporting in solvers

### 3. Cohesive Zone Models (NEW in rc2)

Port of cohesive zone modeling (CZM) from the legacy `czm` branch to `main`:

- `CohesiveZoneModelNIFNSL` base class in the kernel
- `BinaryCohesiveNSL` concrete implementation in `mechanics/src/czm/`
- `CohesiveFrictionContact` one-step nonsmooth problem solver
- Internal variables support in `Interaction` for tracking cohesive state
- Position-level formulation and Bullet contact detection integration

### 4. Parallel & GPU Solvers (NEW in rc2)

Experimental high-performance numerics additions:

- **Parallel PGS/Jacobi solvers** for LCP and 2D/3D friction contact using OpenMP and graph coloring
- **GPU (CUDA) solvers** for friction contact problems with sparse block matrix support
- **PETSc integration** for parallel linear algebra and graph solvers
- Configurable via `WITH_CUDA`, PETSc, and OpenMP CMake options

### 5. Python Bindings: Pybind11 Replaces SWIG

- **pybind11 bindings** for numerics, kernel, mechanics, control, and IO
- **Eigen-NumPy zero-copy mapping** for dense and sparse matrices
- **Improved exception translation** and type safety
- **Boundary conditions** exposed in Python (new in rc2)
- **GPU and parallel solver bindings** added (new in rc2)

### 6. C++ Modernization & Code Quality

- **Namespace encapsulation**: kernel classes moved into `siconos::`
- **C++20 features** used throughout new code
- **Removal of `SP::` typedefs** in favor of explicit `std::shared_ptr<T>`
- **NumericsMatrix API refactoring**: cleaner creation/copy API, storage pointer renaming, robust sparse factorization handling
- **Removal of `exit(EXIT_FAILURE)`** from `NumericsMatrix` error paths
- **Pre-commit hooks** introduced: `clang-format`, `ruff`, `ruff-format`, `shellcheck`

### 7. Build System & Dependencies

- **scikit-build** integration for user-defined plugins (replacing `cppimport`)
- **Eigen3**, **pybind11**, and **nlohmann/json** as key dependencies
- **Guix** packaging improvements: `siconos-devel`, PETSc support, pybind11 3
- **Address Sanitizer (ASan)** support in CI
- **Updated Docker images** (Ubuntu 24.04, Intel oneAPI)

---

## What's New in rc2 Compared to 2026.rc1

### Numerics Hardening

- Refactored and robustified `NumericsMatrix` copy/sparse factorization paths
- Preserved user-provided `problem->M` in IPM solvers with proper Cholesky/LU fallback
- Vectorized error computation and cleaned-up error-scaling options for FC3D/GFC3D
- Added `GFC2D_NSN_AC` with AC/JM formulations, automatic warm start, and configurable Jacobian regularization
- Added GFC reformulation mixed-storage sparse fallback and 2D fixes
- Propagated Fischer-Burmeister non-finite errors instead of aborting
- Standardized verbose output across solvers

### vKernel Maturity

- Registry-backed world alongside `storage::make`
- Dynamic add/remove of items through per-item column lists
- Block diagonal W with stored FEM factorizations
- Disk/fem contact forces output to HDF5
- `vview` FEM mesh visualization with `--view-as-glyphs`
- Compile-time and runtime database naming options (`NONOS_COMPILE_TIME_DATABASE`)
- Split `_nonos` bindings over multiple translation units for faster builds

### Mechanics & IO

- Cohesive zone model integration with position-level formulation
- Renamed `NewtonEuler5DR` → `NewtonEulerVelocityAngularVelocityR` and matching contact relation classes
- FEM strain/stress and constraint tensor output
- `mechanics_run` improvements for `vnative_options`, `contact_nodes`, and empty nodal forces
- Output of cohesive forces and displacement jumps

### Control & Kernel

- Added Lagrangian sparse scleronomous and rheonomous relations (modeling only)
- Fixed shared/weak pointer cycle in simulation graph (ASan)
- Boundary conditions exposed in pybind11 Python API
- Improved `RotationQuaternion` documentation and naming

### Developer Experience

- Added pre-commit configuration with formatting and linting hooks
- Fixed C++ and Python API documentation generation (broken since pybind11 migration)
- Removed most compiler warnings across kernel, numerics, and control

---

## Breaking Changes

### API Changes

1. **Include paths**: Internal numerics headers have been restructured
   - Old: `#include "utils/solver_registry.h"`
   - New: `#include "solver_registry.h"`

2. **Type names** (Plasticity):
   - Old: `MohrCoulomb2DProblem* problem`
   - New: `Plasticity2DProblem* problem`

3. **Enum values** (Plasticity):
   - Old: `SICONOS_MC2D_NSGS`, `MOHR_COULOMB_2D_NSGS`
   - New: `SICONOS_PLASTICITY_2D_NSGS`, `PLASTICITY_2D_NSGS`

4. **Function names** (Plasticity):
   - Old: `mc2d_driver()`, `mohrCoulomb2DProblem_new()`
   - New: `plasticity_2d_driver()`, `plasticity2DProblem_new()`

5. **Relation class rename**:
   - Old: `NewtonEuler5DR`, `Contact5DR`, `Bullet5DR`
   - New: `NewtonEulerVelocityAngularVelocityR`, `ContactVelocityAngularVelocityR`, `BulletVelocityAngularVelocityR`

6. **GFC reformulation functions**:
   - Old: `globalFrictionContactProblem_reformulate_*` / `*_restore`
   - New: `globalFrictionContactProblem_reduce_*` / `expand_*` (deprecated aliases kept)

7. **Numerics public API cleanup**:
   - Expanded cryptic abbreviations in FrictionContact/GlobalFrictionContact solver option setters
     - `fc3d_dsfp_set_default` → `fc3d_de_saxce_fixed_point_set_default`
     - `fc3d_fpp_set_default` → `fc3d_fixed_point_projection_set_default`
     - `fc3d_hp_set_default` → `fc3d_hyperplane_projection_set_default`
     - `fc3d_eg_set_default` → `fc3d_extra_gradient_set_default`
     - `fc3d_pfp_set_default` → `fc3d_panagiotopoulos_fixed_point_set_default`
     - `fc3d_aclmfp_set_default` → `fc3d_alart_curnier_lm_fixed_point_set_default`
     - `fc3d_tfp_set_default` → `fc3d_tresca_fixed_point_set_default`
     - `fc3d_nsn_ac_set_default` → `fc3d_nonsmooth_newton_alart_curnier_set_default`
     - `fc3d_nsn_fb_set_default` → `fc3d_nonsmooth_newton_fischer_burmeister_set_default`
     - `fc3d_nsn_nm_set_default` → `fc3d_nonsmooth_newton_jean_moreau_set_default`
     - `fc3d_nsn_ac_new_set_default` → `fc3d_nonsmooth_newton_alart_curnier_reformulated_set_default`
     - `fc3d_ipm_snm_set_default` → `fc3d_IPM_SNM_set_default`
     - `gfc3d_nsn_ac_set_default` → `gfc3d_nonsmooth_newton_alart_curnier_set_default`
     - `gfc3d_admm_set_default` → `gfc3d_alternating_direction_method_of_multipliers_set_default`
     - `gfc3d_aclmfp_set_default` → `gfc3d_alart_curnier_lm_fixed_point_set_default`
     - `gfc3d_ipm_snm_set_default` → `gfc3d_IPM_SNM_set_default`
   - Renamed misleading `_new` suffixes
     - `computeAlartCurnier2D_new` → `computeAlartCurnier2D_reformulated`
     - `fc3d_nonsmooth_Newton_AlartCurnier_new` → `fc3d_nonsmooth_Newton_AlartCurnier_reformulated`
     - `gfc3d_nonsmooth_Newton_AlartCurnier_new_wr` → `gfc3d_nonsmooth_Newton_AlartCurnier_reformulated_wr`
     - `NM_internalData_new` → `NM_internalData_alloc`
   - Standardized problem struct lifecycle verbs
     - `freeAVI` → `AVI_free`
     - `freeNCP` → `NCP_free`
   - Deprecated aliases are centralized in `numerics/src/utils/deprecated_numerics_api.h`

8. **Numerics public API cleanup — wave 2**:
   - Expanded remaining cryptic abbreviations
     - `gfc3d_vi_eg_set_default` → `gfc3d_VI_extra_gradient_set_default`
     - `gfc3d_vi_fpp_set_default` → `gfc3d_VI_fixed_point_projection_set_default`
     - `gfc2d_vi_eg_set_default` → `gfc2d_VI_extra_gradient_set_default`
     - `fc3d_FixedP_initialize` → `fc3d_fixed_point_initialize`
     - `fc3d_FixedP_solve` → `fc3d_fixed_point_solve`
     - `fc3d_FixedP_free` → `fc3d_fixed_point_free`
     - `fc3d_FixedP_computeError` → `fc3d_fixed_point_compute_error`
     - `F_GlockerFixedP` → `F_Glocker_fixed_point`
     - `compute_Z_GlockerFixedP` → `compute_Z_Glocker_fixed_point`
   - Standardized `FrictionContactProblem` and `GlobalFrictionContactProblem` lifecycle/display verbs
     - `frictionContact_display` → `frictionContactProblem_display`
     - `frictionContact_printInFile` → `frictionContactProblem_printInFile`
     - `frictionContact_printInFilename` → `frictionContactProblem_printInFilename`
     - `globalFrictionContact_display` → `globalFrictionContactProblem_display`
     - `globalFrictionContact_free` → `globalFrictionContactProblem_free`
     - `globalFrictionContact_printInFile` → `globalFrictionContactProblem_printInFile`
     - `globalFrictionContact_printInFileName` → `globalFrictionContactProblem_printInFilename`
   - Removed dead NCP error helpers
     - `Compute_NCP_error1`
     - `Compute_NCP_error2`

### Build System

- **SWIG is being phased out**; pybind11 is the future for Python bindings
- **New required dependencies**:
  - Eigen3 (for vKernel module)
  - pybind11 (for Python bindings)
  - Boost headers (for utilities, not ublas)
  - nlohmann/json
- **Optional new dependencies**:
  - CUDA toolkit (for GPU solvers)
  - PETSc (for parallel solvers)
- **Scikit-build** for plugin development (replacing `cppimport`)

---

## Migration Guide

### For C API Users

1. Update include paths for internal numerics headers (if using them)
2. Replace `MohrCoulomb2DProblem` with `Plasticity2DProblem`
3. Update solver enum names if using plasticity solvers
4. Update GFC reformulation function names if used directly
5. Update renamed numerics function names or include `numerics/src/utils/deprecated_numerics_api.h` for backward-compatible aliases

### For C++ Users

1. Add `using namespace siconos;` or qualify types with `siconos::`
2. Replace `SP::Type` with `std::shared_ptr<Type>`
3. Update any custom solver integration to use the registration system
4. Rename `NewtonEuler5DR` and related contact relations to `NewtonEulerVelocityAngularVelocityR`

### For Python Users

1. No immediate changes required for legacy SWIG-based bindings
2. pybind11-based bindings are the future; minor API adjustments may be needed when migrating

---

## Testing

- Core numerics and kernel tests passing in standard, ASan, and hard-debug/sanitize modes
- New test collections for GFC2D, rolling friction, cohesive friction, and vKernel FEM/disk scenarios
- GPU solver tests added (conditional on CUDA)
- Parallel PGS/Jacobi tests added (conditional on OpenMP/PETSc)
- Pre-commit hooks applied across Python, shell, and C++ code

---

## Contributors

This release includes contributions from the Siconos development team at INRIA and external collaborators.

---

## Full Changelog

See the [Changelog](../../Changelog) file for a detailed list of changes.

For the previous release candidate notes, see [RELEASE_2026.rc1.md](RELEASE_2026.rc1.md).

---

## Download & Installation

```bash
# Clone the repository
git clone https://gricad-gitlab.univ-grenoble-alpes.fr/nonsmooth/siconos.git
cd siconos
git checkout 2026.rc2

# Build with CMake
mkdir build && cd build
cmake .. -DCMAKE_BUILD_TYPE=Release
make -j$(nproc)
```

Or use Docker:
```bash
docker run --rm -p 8888:8888 -ti \
  gricad-registry.univ-grenoble-alpes.fr/nonsmooth/siconos-tutorials/siconoslab-main:latest
```

---

## Documentation & Support

- **Documentation**: https://nonsmooth.gricad-pages.univ-grenoble-alpes.fr/siconos
- **Tutorials**: https://nonsmooth.gricad-pages.univ-grenoble-alpes.fr/siconos-tutorials
- **Issue Tracker**: https://gricad-gitlab.univ-grenoble-alpes.fr/nonsmooth/siconos/-/issues
- **Homepage**: http://siconos.org

---

**License**: Apache License 2.0
**Copyright**: 2024-2026 INRIA
