# Siconos 2026.rc1 Release Notes

**Release Date:** February 2026
**Previous Release:** 4.5.0
**Commits:** 1193+ changes across the entire codebase

---

## Overview

Siconos 2026.rc1 is a major release featuring a **brand new vKernel architecture** for large-scale systems, significant modernization of the C++ codebase, a complete refactoring of the numerics solver infrastructure, and new Python bindings. This release represents a substantial evolution of the software architecture while maintaining backward compatibility where possible.

---

## Major Highlights

### 1. vKernel: New Architecture for Large-Scale Systems (BRAND NEW)

Introduction of **vKernel** (virtual kernel), a completely new computational backend designed for high-performance simulation of large-scale nonsmooth systems:

#### Structure of Arrays (SoA) Architecture

The cornerstone of vKernel is its **SoA (Structure of Arrays)** data layout, replacing the traditional AoS (Array of Structures) approach:

- **Cache-friendly memory layout**: Data stored as contiguous arrays of scalars rather than arrays of objects
  - Example: All positions `[x1, y1, z1, x2, y2, z2, ...]` stored contiguously instead of `[Object1(x,y,z), Object2(x,y,z), ...]`
- **Vectorized operations**: Enables SIMD-friendly access patterns for modern CPUs
- **Better branch prediction**: Homogeneous data types allow for more predictable execution
- **Reduced memory footprint**: Eliminates padding and pointer overhead from object-oriented structures

#### Key Components

- **SoA Storage System** (`vkernel/soa/`):
  - `storage/` - Generic SoA container implementations
  - `topology/` - Graph-based system topology using SoA patterns
  - `algebra/` - Linear algebra operations optimized for SoA layout

- **Type-Safe Algebra** (`vkernel/soa/src/siconos/algebra/eigen.hpp`):
  - Template-based matrix/vector concepts with compile-time size checking
  - Unified interface for sparse and dense operations
  - Zero-copy memory mapping between C++ and Python

- **Runtime-Defined Systems**:
  - Dynamical systems and interactions defined at runtime rather than compile-time
  - Flexible topology management using modern C++20 features
  - Attached storage for user-defined properties

#### Status
The vKernel module is available as an **experimental parallel architecture** to the classic kernel. It represents the future direction of Siconos for large-scale simulations (millions of objects). Future releases will continue expanding vKernel capabilities and eventually make it the primary backend.

**Key directories:**
- `vkernel/soa/src/siconos/algebra/` - Eigen-based algebra
- `vkernel/soa/src/siconos/storage/` - SoA containers
- `vkernel/soa/src/siconos/topology/` - Graph topology management

---

### 2. Solver Registration System & Unified Drivers (Numerics)

A complete redesign of the solver infrastructure in the numerics module:

- **New solver registration system** (`solver_registry.h`): Solvers now self-register using `REGISTER_SOLVER` and `REGISTER_SOLVER_3VAR` macros, enabling runtime solver lookup and dispatch
- **Unified driver architecture**: All drivers now use the solver registry for dispatch, eliminating large switch statements
- **Standardized function signatures**: Consistent `init/solve/free` patterns across all solver families
- **Support for 3-variable solvers**: Proper handling of GFC3D and GRFC3D solvers that require `globalVelocity`

**Affected solver families:**
- FrictionContact (fc3d, fc2d)
- GlobalFrictionContact (gfc3d, gfc2d)
- RollingFrictionContact (new solvers added, files renamed to `rolling_friction_3d_*`)
- GlobalRollingFrictionContact (new solvers added, files renamed to `global_rolling_friction_3d_*`)
- Plasticity (files renamed from `mc2d_*` to `plasticity_2d_*`, types renamed from `MohrCoulomb2DProblem` to `Plasticity2DProblem`)
- LCP, MLCP, MCP, AVI, VI, Relay, SOC-LCP

---

### 3. Python Bindings: Pybind11 Replaces SWIG

Complete overhaul of Python bindings architecture:

- **Pybind11 replaces SWIG**: Modern C++11-based binding generation with better type safety and performance
- **Eigen-Python mapping**: Direct mapping between Eigen matrices and NumPy arrays (zero-copy when possible)
- **Improved error messages**: Better exception translation between C++ and Python
- **Simplified build**: No more SWIG interface files (`.i`), bindings are pure C++

**Modules with pybind11 bindings:**
- **Numerics**: FrictionContact, LCP, MLCP, Relay, MCP, AVI
- **Kernel**: Core classes, algebra, modeling, simulation
- **Mechanics**: Contact detection, relations, joints
- **Control**: Sliding mode controllers, observers
- **IO**: HDF5 serialization

**Key features:**
- `py::class_` and `py::classh` for memory management
- `eigen2python_pybind11` for matrix/vector conversions
- Support for both dense and sparse matrix transfers

---

### 4. C++ Modernization & Architecture

Major C++ codebase modernization:

- **Namespace encapsulation**: All kernel classes moved into `siconos::` namespace
- **C++20 features**: Concepts, structured bindings, range-based algorithms where applicable
- **Removal of SP:: typedefs**: Using explicit `std::shared_ptr<T>` instead
- **Constructor delegation**: Cleaner initialization patterns
- **Enum class**: Type-safe enumerations replacing plain enums
- **Template improvements**: Better generic programming support

---

### 5. Build System & Dependencies

- **Scikit-build integration**: Modern Python packaging for user-defined functions (plugins), replacing cppimport
- **Eigen3**: New dependency for vKernel linear algebra (replaces Boost::ublas in new code)
- **pybind11**: New dependency for Python bindings (replacing SWIG)
- **Boost**: Now used only for utilities, not for ublas
- **nlohmann/json**: Integrated for configuration and data exchange
- **CMake improvements**: Better component management, stricter warning levels
- **CI/CD updates**: New Docker images (Ubuntu 24.04, Intel oneAPI), updated GitLab CI

---

### 6. New Utilities & Infrastructure

- **Tolerance Manager** (`tolerance_manager.h`): Unified tolerance adaptation across iterative solvers
- **Error Computation Utilities** (`error_computation.h`): Standardized error metrics
- **Safe Integer Conversions** (`safe_int_conversion.h`): Type-safe conversions between BLAS/LAPACK integers and C integers
- **Matrix Versioning**: Optimized matrix update tracking

---

### 7. Bug Fixes & Quality Improvements


- **Memory safety**: Fixed various memory leaks and initialization issues
- **Include structure flattened**: Internal headers now use flat paths (e.g., `"solver_registry.h"` instead of `"utils/solver_registry.h"`)
- **Warning cleanup**: Removed pragma warnings, fixed strict prototype warnings

---

## New Features Summary

| Category | Feature | Status |
|----------|---------|--------|
| **vKernel** | SoA (Structure of Arrays) architecture | 🔄 Experimental |
| **vKernel** | Runtime-defined dynamical systems | 🔄 Experimental |
| **vKernel** | Eigen linear algebra backend | 🔄 Experimental |
| **Numerics** | Solver registration system | ✅ Complete |
| **Numerics** | Rolling friction 3D solvers | ✅ Complete |
| **Numerics** | Global rolling friction 3D (IPM, NSGS-WR) | ✅ Complete |
| **Numerics** | Plasticity 2D (renamed from Mohr-Coulomb) | ✅ Complete |
| **Numerics** | Unified NSGS framework | ✅ Complete |
| **Numerics** | Tolerance manager | ✅ Complete |
| **Kernel** | Namespace encapsulation | ✅ Complete |
| **Kernel** | C++20 migration | 🔄 In Progress |
| **Python** | pybind11 bindings (replacing SWIG) | 🔄 In Progress |
| **Python** | Eigen-NumPy zero-copy mapping | ✅ Complete |
| **Build** | scikit-build support | ✅ Complete |

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

### Build System

- **SWIG is being phased out**; pybind11 is the future for Python bindings
- **New required dependencies**:
  - Eigen3 (for vKernel module)
  - pybind11 (for Python bindings)
  - Boost headers (for utilities, not ublas)
  - nlohmann/json
- Updated CMake minimum requirements for some features
- **Scikit-build** for plugin development (replacing cppimport)

---

## Migration Guide

### For C API Users

1. Update include paths for internal numerics headers (if using them)
2. Replace `MohrCoulomb2DProblem` with `Plasticity2DProblem`
3. Update solver enum names if using plasticity solvers

### For C++ Users

1. Add `using namespace siconos;` or qualify types with `siconos::`
2. Replace `SP::Type` with `std::shared_ptr<Type>`
3. Update any custom solver integration to use the registration system

### For Python Users

1. No immediate changes required for SWIG-based bindings
2. Future releases will transition to pybind11 bindings

---

## Testing

- **96 numerics tests** passing (100% of core numerics test suite)
- New test collections for rolling friction contact
- Updated plasticity test suite with renamed data files
- ASan (Address Sanitizer) support added to CI

---

## Contributors

This release includes contributions from the Siconos development team at INRIA and external collaborators.

---

## Full Changelog

See the [Changelog](./Changelog) file for a detailed list of changes.

---

## Download & Installation

```bash
# Clone the repository
git clone https://gricad-gitlab.univ-grenoble-alpes.fr/nonsmooth/siconos.git
cd siconos
git checkout 2026.rc1

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
