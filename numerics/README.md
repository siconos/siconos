# Siconos Numerics

This directory contains the **Siconos Numerics** library - a collection of numerical algorithms for solving non-smooth dynamical systems and complementarity problems.

## Overview

Siconos Numerics provides solvers for various problem formulations including:

- **Friction Contact Problems** - Coulomb friction in 2D and 3D
- **Rolling Friction Contact Problems** - Friction with rolling resistance
- **Complementarity Problems** - LCP, NCP, MLCP, MCP
- **Variational Inequalities** - AVI, VI
- **Optimization Problems** - QP, SOCP

## Directory Structure

| Directory | Purpose |
|-----------|---------|
| `src/` | Source code for all numerical algorithms |
| `python/` | Python bindings (pybind11) |
| `share/` | Shared data files |

## Detailed Documentation

**For detailed information about the source code organization, solver naming conventions, and migration guides, see:**

### [src/README.md](src/README.md)

The `src/README.md` file contains:
- Complete directory structure description
- Problem formulation directories (FrictionContact, RollingFrictionContact, LCP, etc.)
- Generic solver infrastructure (NonSmoothSolvers)
- Solver naming conventions and ID ranges
- Key architectural patterns
- Migration guides for recent refactoring changes

## Quick Links

- [Build Instructions](../INSTALL) - How to build Siconos Numerics
- [Source Documentation](src/README.md) - Detailed source organization
- [Tests](numerics_tests.cmake) - Test configuration
- [Doxygen Documentation](../../docs/) - API reference

## License

See [COPYING](../COPYING) for license information.
