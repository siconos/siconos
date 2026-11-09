# Siconos Kernel

This directory contains the **Siconos Kernel** library — the core C++ simulation engine for modeling and integrating nonsmooth dynamical systems.

## Overview

The Kernel provides the high-level abstractions used to describe a nonsmooth dynamical system and to simulate it over time. It bridges the modeling layer (dynamical systems, interactions, relations, and nonsmooth laws) with the numerical solvers from [`numerics`](../numerics).

Key responsibilities of the Kernel include:

- **Algebra**: vectors, matrices, and block objects used throughout Siconos.
- **Modeling**: dynamical systems, relations, interactions, and nonsmooth laws.
- **Simulation**: time discretization, one-step integrators, one-step nonsmooth problems, and the main simulation loop.
- **Utilities**: visitors, serialization helpers, exception handling, and graph/property utilities.

## Directory Structure

| Directory | Purpose |
|-----------|---------|
| `src/` | Main C++ source code |
| `src/algebra/` | Linear algebra: `SiconosVector`, `SiconosMatrix`, `BlockVector`, `SiconosMemory`, etc. |
| `src/modeling/` | Dynamical systems, relations, interactions, and nonsmooth laws |
| `src/simulation/` | Simulation classes, one-step integrators (OSI), one-step nonsmooth problems (OSNS), events manager |
| `src/utils/` | Common utilities: exceptions, graph tools, serialization, quaternion helpers |
| `python/` | Python bindings using pybind11 |
| `swig/` | SWIG-related test support (legacy) |
| `tests-common/` | Shared testing utilities for the Kernel |

## Key Modeling Components

### Dynamical Systems

- `DynamicalSystem` — base class
- `FirstOrderLinearDS`, `FirstOrderNonLinearDS` — first-order systems
- `LagrangianDS`, `LagrangianLinearTIDS`, `LagrangianSparseDS` — Lagrangian (second-order) systems
- `NewtonEulerDS` — 3D Newton-Euler rigid-body systems

### Relations and Interactions

- `Relation` — base class for relations between variables
- `Interaction` — links a relation, a nonsmooth law, and dynamical systems
- Lagrangian and first-order relation families (`LagrangianScleronomousR`, `FirstOrderType1R`, ...)

### Nonsmooth Laws

- `NewtonImpactNSL`, `NewtonImpactFrictionNSL`, `NewtonImpactRollingFrictionNSL`
- `RelayNSL`, `NormalConeNSL`, `MultipleImpactNSL`
- `ComplementarityConditionNSL`, `MixedComplementarityConditionNSL`

## Simulation Components

### Simulation Classes

- `Simulation` — abstract base
- `TimeStepping`, `TimeSteppingCombinedProjection`, `TimeSteppingDirectProjection`, `TimeSteppingD1Minus` — time-stepping schemes
- `EventDriven` — event-driven schemes

### One-Step Integrators (OSI)

- `MoreauJeanOSI`, `MoreauJeanGOSI`, `MoreauJeanBilbaoOSI`, `MoreauJeanCombinedProjectionOSI`, `MoreauJeanDirectProjectionOSI`
- `EulerMoreauOSI`, `NewMarkAlphaOSI`, `ZeroOrderHoldOSI`
- `LsodarOSI`, `Hem5OSI`, `D1MinusLinearOSI`, `SchatzmanPaoliOSI`

### One-Step Nonsmooth Problems (OSNS)

- `LCP`, `MLCP`, `MLCPProjectOnConstraints`
- `FrictionContact`, `RollingFrictionContact`, `GlobalFrictionContact`, `GlobalRollingFrictionContact`, `CohesiveFrictionContact`
- `AVI`, `Relay`, `QP`, `Equality`, `MultipleImpact`, `GenericMechanical`

## Python Bindings

Python bindings are generated with pybind11 in the `python/` directory:

- `modeling/` — dynamical systems, relations, nonsmooth laws
- `simulation/` — simulation, OSI, OSNS
- `integrators/` — integrator bindings
- `nonsmooth_formulations/` — nonsmooth problem formulations
- `geometry_pybind11.cpp` — geometry helpers

## Tests

Unit tests are located alongside the source code in `src/<module>/test/` and additional Python tests live in `python/tests/`. The test suite is configured in [`kernel_tests.cmake`](kernel_tests.cmake).

## Dependencies

The Kernel links against:

- [`externals`](../externals) — bundled external libraries
- [`numerics`](../numerics) — nonsmooth problem solvers
- Boost
- Eigen3
- GMP
- Optional: nlohmann JSON (`WITH_JSON`), Boost.Log (`WITH_BOOST_LOG`)

## Documentation

- [Siconos documentation](https://nonsmooth.gricad-pages.univ-grenoble-alpes.fr/siconos)
- [Architecture overview](../docs/sphinx/devel_guide/architecture.rst)
- [Build instructions](../INSTALL)
- [Main README](../README.md)

## License

See [COPYING](../COPYING) for license information.
