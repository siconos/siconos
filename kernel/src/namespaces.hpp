// SPDX-FileCopyrightText: 2026 INRIA
//
// SPDX-License-Identifier: Apache-2.0
//
// Siconos is a program dedicated to modeling, simulation and control
// of non smooth dynamical systems.

/** @file namespaces.hpp description of all siconos namespaces
 * Only useful for doc.
 */

/** @namespace siconos
 * @brief Main siconos namespace
 */
namespace siconos {}

/** @namespace siconos::exception
 *  @brief Exceptions handling
 */
namespace siconos::exception {}

/** @namespace siconos::algebra
 *  @brief Classes and tools to handle matrices and vectors definitions and operations in
 *  siconos
 */
namespace siconos::algebra {

/** @namespace siconos::algebra::io
 * @brief utilities to handle file input/output for vectors and matrices
 */
namespace io {}

}  // namespace siconos::algebra

/** @namespace siconos::modeling
 *  @brief Classes and tools used to define and configure the nonsmooth dynamical system(s)
 *
 */
namespace siconos::modeling {}

/** @namespace siconos::graphs
 *  @brief Definition and description of the graph structures (ds and interactions graphs
 * mostly)
 */
namespace siconos::integrators {}

/** @namespace siconos::integrators
 *  @brief Classes and tools for time integration of the dynamics (one-step integrators)
 */
namespace siconos::integrators {}

/** @namespace siconos::nonsmooth_formulations
 *  @brief Classes and tools dedicated to the description of the nonsmooth problem (formulation
 * and solver)
 */
namespace siconos::nonsmooth_formulations {}

/** @namespace siconos::simulation
 *  @brief Classes and tools dedicated to the description of the simulation of a NSDS (time
 * discretization, ...)
 */
namespace siconos::simulation {}
