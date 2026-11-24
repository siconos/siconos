/* Siconos is a program dedicated to modeling, simulation and control
 * of non smooth dynamical systems.
 *
 * SPDX-FileCopyrightText: INRIA
 * SPDX-License-Identifier: Apache-2.0
 */

/*! \file SiconosConst.hpp
\brief General constants for Siconos Kernel.
*/

#ifndef KernelConst_H
#define KernelConst_H

#include <limits>

// const or tools for internal use (i.e. not intended to be in the user API)
namespace siconos::internal {

/**
 Internal bound max levels for time integrators.
 This value may be checked to see if initialization has occured.
*/
constexpr auto LEVELMAX = 999;

/** double precision machine */
constexpr double MACHINE_PREC = std::numeric_limits<double>::epsilon();

/** default tolerance for simulation algorithms. */
constexpr double DEFAULT_TOLERANCE = 10 * MACHINE_PREC;

/** default tolerance for EventDriven algorithms */
constexpr double DEFAULT_TOL_ED = 1000 * DEFAULT_TOLERANCE;

/** tick default value (for events in event-driven scheme)
 *  it has to be greater than DBL_EPSILON */
constexpr double DEFAULT_TICK = 1e-16;

// Events management stuff
constexpr unsigned long int GAPLIMIT_DEFAULT = 100;

}  // namespace siconos::internal

#endif
