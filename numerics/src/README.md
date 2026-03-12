# Siconos Numerics Source Directory

This directory contains the core numerical algorithms of the Siconos library for solving non-smooth dynamical systems.

## Directory Structure

### Problem Formulation Directories

These directories contain problem-specific data structures, solvers, and utilities organized by mathematical formulation.

#### `FrictionContact/` - Friction Contact Problems
Contact mechanics with Coulomb friction in 2D and 3D.
- **Problem types**: 2D/3D frictional contact, global formulations
- **Key solvers**: NSGS, Non-Smooth Newton (Alart-Curnier, Fischer-Burmeister), Proximal, ADMM, IPM
- **Files**: `FrictionContact_options.h` (solver IDs and options)

#### `LCP/` - Linear Complementarity Problems
Linear complementarity between a variable and its dual.
- **Solvers**: Lemke, Pivoting methods, Newton-based methods

#### `MCP/` - Mixed Complementarity Problems
Mixed complementarity combining equality and inequality constraints.

#### `MLCP/` - Mixed Linear Complementarity Problems
Linear complementarity with additional equality constraints.

#### `NCP/` - Nonlinear Complementarity Problems
Nonlinear variants of complementarity problems.

#### `AVI/` - Affine Variational Inequalities
Affine variational inequality problems.

#### `VI/` - Variational Inequalities
General variational inequality formulations.
- **Solvers**: Extra-gradient, Fixed-point projection, Hyperplane projection

#### `Relay/` - Relay Systems
Systems with relay (switching) non-smoothness.

#### `RollingFrictionContact/` - Rolling Friction Contact Problems
Contact mechanics with Coulomb friction and rolling resistance in 2D and 3D.
- **Problem types**: 2D/3D rolling friction contact, global rolling formulations
- **Key solvers**: NSGS, ADMM, IPM
- **Files**: 
  - `RollingFrictionContact_options.h` (solver IDs and options)
  - `RollingFrictionContactProblem.h/.c` - Local rolling contact problem
  - `GlobalRollingFrictionContactProblem.h/.c` - Global rolling contact problem

**Note**: Rolling friction was previously in `FrictionContact/` but has been moved to its own directory for better organization.

#### `SOCP/` - Second-Order Cone Programming
Conic optimization problems.

#### `QP/` - Quadratic Programming
Convex quadratic programming problems.

### Generic Solver Infrastructure

#### `NonSmoothSolvers/` - Generic Non-Smooth Solvers
Framework for problem-agnostic non-smooth solvers.

| File | Purpose |
|------|---------|
| `NonSmoothGaussSeidel_options.h` | NSGS algorithm options (shuffle, freezing, relaxation, error evaluation) |
| `nsgs_generic.h` | Generic NSGS framework (single-loop and 4-stage variants) |
| `nsgs_generic_driver.h` | Driver utilities for NSGS solvers |
| `nsgs_generic_instrumented.h` | Profiled/instrumented NSGS variants |
| `NonSmoothNewton.h/.c` | Generic Non-Smooth Newton methods |
| `Newton_methods.h/.c` | Line-search algorithms (Armijo, Goldstein) for Newton |

**Design Principle**: These solvers work with any problem type through callback functions. The NSGS framework is used by:
- Friction contact problems (2D, 3D)
- Rolling friction contact problems (2D, 3D)
- Mohr-Coulomb plasticity problems
- Generic mechanical systems

### Supporting Infrastructure

#### `core/` - Core Types and Options
**NEW**: Fundamental types used across all modules.
- `SolverOptions.h/.c` - Solver configuration structure
- `SiconosSets.h/.c` - Set operations
- `NumericsDataVersion.h/.c` - Version tracking
- `sn_error_handling.h/.c` - Error handling

#### `linalg/` - Linear Algebra
**NEW**: Matrix operations and linear algebra utilities.
- `NumericsMatrix.h/.c` - Matrix representations (dense, sparse, block)
- `NumericsVector.h/.c` - Vector operations
- `NumericsSparseMatrix.h/.c` - Sparse matrix formats
- `CSparseMatrix.h/.c` - CSparse integration
- `SparseBlockMatrix.h/.c` - Block-structured matrices
- `SiconosBlas.h/.c` - BLAS wrappers
- `SiconosLapack.h/.c` - LAPACK wrappers
- `op3x3.h/.c`, `op5x5.h` - Small matrix optimizations
- `projectionOnCone.h/.c`, etc. - Projections onto cones

#### `optimization/` - Optimization Methods
**NEW**: Line search and optimization utilities.
- `line_search.h/.c` - Line search framework
- `ArmijoSearch.h/.c` - Armijo backtracking
- `GoldsteinSearch.h/.c` - Goldstein line search
- `NMS.h/.c` - Non-monotone line search
- `FischerBurmeister.h/.c` - Fischer-Burmeister functions
- `JordanAlgebra.h/.c` - Jordan algebra operations
- `min_merit.h/.c`, `Qi_merit.h/.c` - Merit functions

#### `io/` - Input/Output
**NEW**: File I/O and logging utilities.
- `hdf5_logger.h/.c` - HDF5 logging
- `io_tools.h/.c` - General I/O utilities
- `sn_logger.h` - Siconos numerics logging
- `GAMSlink.h/.c` - GAMS interface

#### `utils/` - General Utilities
**NEW**: Miscellaneous utility functions.
- `siconos_debug.h/.c` - Debug utilities
- `NumericsVerbose.h` - Verbose output
- `enum_tool.h/.c` - Enum utilities
- `graph.h/.c` - Graph algorithms
- `quartic.h/.c` - Quartic equation solver
- `vertex_extraction.h/.c` - Vertex extraction

#### `tools/` - Backward Compatibility Layer
**LEGACY**: This directory now serves as a backward compatibility layer. Headers here are wrappers that include files from their new locations (`core/`, `linalg/`, `optimization/`, `io/`, `utils/`).

**Note**: Please update your includes to use the new paths. The wrapper headers may be deprecated in a future release.

#### `GenericMechanical/` - Multi-Physics Systems
Coupled systems combining multiple problem types (friction + impact + ...).

#### `Plasticity/` - Plasticity Models
Mohr-Coulomb and other plasticity models for granular materials.

### Deprecated

#### `Unstable_or_deprecated/`
Legacy code maintained for backward compatibility but not recommended for new development.

## Solver Naming Conventions

### Solver IDs
Solver identifiers follow the pattern:
```
SICONOS_<FORMULATION>_<METHOD>[_WR]
```

Examples:
- `SICONOS_FRICTION_3D_NSGS` - Non-Smooth Gauss-Seidel for 3D friction contact
- `SICONOS_FRICTION_3D_NSN_AC` - Non-Smooth Newton (Alart-Curnier) for 3D friction
- `SICONOS_GLOBAL_FRICTION_3D_ADMM_WR` - ADMM for global formulation with wrapper
- `SICONOS_ROLLING_FRICTION_3D_NSGS` - NSGS for 3D rolling friction contact
- `SICONOS_GLOBAL_ROLLING_FRICTION_3D_IPM` - IPM for global rolling friction contact

### Parameter Indices
Configuration options use typed indices:

| Prefix | Array | Purpose |
|--------|-------|---------|
| `SICONOS_<FORMULATION>_IPARAM_*` | `iparam[]` | Integer options (max iter, strategies) |
| `SICONOS_<FORMULATION>_DPARAM_*` | `dparam[]` | Double options (tolerances, parameters) |

## Key Architectural Patterns

### 1. Problem-Specific + Generic Solvers
Each problem formulation implements its own solvers AND can use generic frameworks:
- Direct solvers: `fc3d_nsgs.c` (friction-specific NSGS)
- Generic framework: Uses `nsgs_generic.h` with callbacks

### 2. Local Problem Solvers
Many solvers (especially NSGS) use internal "one-contact" solvers:
- `SICONOS_ONECONE_NSN` - Newton for single contact
- `SICONOS_ONECONE_ProjectionOnCone` - Projection for single contact

### 3. Solver Composition
Global solvers can wrap local solvers:
- `SICONOS_GLOBAL_FRICTION_3D_NSGS_WR` - NSGS on global formulation

## Migration Guide

### Recent Changes (Refactoring)

#### 1. NSGS Options Moved
Options for shuffle, freezing, relaxation are now in `NonSmoothSolvers/NonSmoothGaussSeidel_options.h`

**Old include:**
```c
#include "FrictionContact_options.h"
// Uses SICONOS_FRICTION_3D_NSGS_SHUFFLE, etc.
```

**New include (optional):**
```c
#include "NonSmoothSolvers/NonSmoothGaussSeidel_options.h"
// Uses SICONOS_NSGS_SHUFFLE, etc.
```

**Backward compatibility:** Old names (`SICONOS_FRICTION_3D_NSGS_*`) remain as macros.

#### 2. Tools Directory Split
Files from `tools/` have been moved to specialized directories:

| Old Location | New Location |
|--------------|--------------|
| `tools/SolverOptions.h` | `core/SolverOptions.h` |
| `tools/NumericsMatrix.h` | `linalg/NumericsMatrix.h` |
| `tools/line_search.h` | `optimization/line_search.h` |
| `tools/hdf5_logger.h` | `io/hdf5_logger.h` |
| `tools/siconos_debug.h` | `utils/siconos_debug.h` |

**Backward compatibility:** Wrapper headers in `tools/` maintain old includes.

#### 3. Enum Naming Standardized
Enums no longer use redundant `_ENUM` suffix:

**Old:**
```c
enum SICONOS_FRICTION_3D_ADMM_IPARAM_ENUM { ... };
```

**New:**
```c
enum SICONOS_FRICTION_3D_ADMM_IPARAM { ... };
```

**Backward compatibility:** Macros map old names to new names.

#### 4. Rolling Friction Contact Moved
Rolling friction contact code has been moved from `FrictionContact/` to its own directory:

**Old includes:**
```c
#include "FrictionContact/RollingFrictionContactProblem.h"
#include "FrictionContact/GlobalRollingFrictionContactProblem.h"
#include "FrictionContact/rolling_fc_Solvers.h"
```

**New includes:**
```c
#include "RollingFrictionContact/RollingFrictionContactProblem.h"
#include "RollingFrictionContact/GlobalRollingFrictionContactProblem.h"
#include "RollingFrictionContact/rolling_fc_Solvers.h"
```

**Solver IDs**: Rolling friction solver IDs have moved to `RollingFrictionContact_options.h`:
- `SICONOS_ROLLING_FRICTION_3D_NSGS = 3000`
- `SICONOS_ROLLING_FRICTION_3D_ADMM = 3003`
- `SICONOS_ROLLING_FRICTION_2D_NSGS = 4000`
- `SICONOS_GLOBAL_ROLLING_FRICTION_3D_NSGS_WR = 5000`
- `SICONOS_GLOBAL_ROLLING_FRICTION_3D_IPM = 5001`

**Backward compatibility:** Wrapper headers in `FrictionContact/` maintain old includes with deprecation warnings.

### For New Solvers
When adding a new solver:
1. Define solver ID in the appropriate problem header
2. For NSGS variants: Use the generic framework in `NonSmoothSolvers/`
3. For Newton variants: Extend `NonSmoothSolvers/Newton_methods.h`
4. Add string name in corresponding `.c` file

## See Also
- `../test/` - Test suites for each problem type
- `../python/` - Python bindings (pybind11)
- `../../docs/` - Doxygen documentation
