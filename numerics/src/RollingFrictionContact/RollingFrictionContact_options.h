// SPDX-FileCopyrightText: INRIA
//
// SPDX-License-Identifier: Apache-2.0
//
// Siconos is a program dedicated to modeling, simulation and control
// of non smooth dynamical systems.

/*!\file RollingFrictionContact_options.h
 * \brief Solver options and constants for Rolling Friction Contact problems
 *
 * This header defines:
 * - Solver IDs (enum ROLLING_FRICTION_SOLVER): All available solvers for rolling
 *   friction contact problems
 * - String names for rolling friction solvers
 *
 * Rolling friction extends standard friction contact with additional rolling
 * resistance constraints. This problem is important in applications involving
 * rolling objects like wheels, balls, or cylindrical contacts.
 *
 * The rolling friction contact solvers are organized into:
 * - 2D Rolling Friction solvers
 * - 3D Rolling Friction solvers (local formulation)
 * - Global formulation solvers
 *
 * Common solver parameters (iparam/dparam) are inherited from FrictionContact_options.h
 */

#ifndef ROLLING_FRICTION_CONTACT_OPTIONS_H
#define ROLLING_FRICTION_CONTACT_OPTIONS_H

#include "FrictionContact_options.h"
#include "NumericsFwd.h"
#include "SolverOptions.h"

/* ===========================================================================
 * Rolling Friction Solver IDs
 * ===========================================================================
 * These extend the FRICTION_SOLVER enum from FrictionContact_options.h
 * with rolling friction specific solvers.
 */

/** \enum ROLLING_FRICTION_SOLVER
 * \brief Unique identifiers for rolling friction contact solvers
 *
 * Solver IDs are organized by problem formulation (2D, 3D, global).
 * IDs 3000-3999: 3D Rolling Friction (local formulation)
 * IDs 4000-4999: 2D Rolling Friction
 * IDs 5000-5999: Global Rolling Friction formulation
 */
enum ROLLING_FRICTION_SOLVER {
  /* -----------------------------------------------------------------------
   * 3D Rolling Friction solvers - Local formulation (IDs 3000-3999)
   * -----------------------------------------------------------------------
   * These solvers work on the local (contact-level) formulation of the
   * 3D rolling friction contact problem.
   */
  /** Non-smooth Gauss-Seidel for 3D rolling friction, local formulation */
  SICONOS_ROLLING_FRICTION_3D_NSGS = 3000,
  /** Projection on cone for 3D rolling friction, one contact solver */
  SICONOS_ROLLING_FRICTION_3D_ONECONTACT_ProjectionOnCone = 3001,
  /** Projection on cone with local iteration for 3D rolling friction */
  SICONOS_ROLLING_FRICTION_3D_ONECONTACT_ProjectionOnConeWithLocalIteration = 3002,
  SICONOS_ROLLING_FRICTION_3D_ONECONTACT_NSN = 3003,
  SICONOS_ROLLING_FRICTION_3D_ONECONTACT_NSN_GP = 3004,
  SICONOS_ROLLING_FRICTION_3D_ONECONTACT_NSN_GP_HYBRID = 3005,

  /** ADMM for 3D rolling friction, local formulation */
  SICONOS_ROLLING_FRICTION_3D_ADMM = 3010,

  /* -----------------------------------------------------------------------
   * 2D Rolling Friction solvers (IDs 4000-4999)
   * -----------------------------------------------------------------------
   * These solvers handle 2D rolling friction contact problems.
   */
  /** Non-smooth Gauss-Seidel for 2D rolling friction */
  SICONOS_ROLLING_FRICTION_2D_NSGS = 4000,
  /** Projection on cone for 2D rolling friction, one contact solver */
  SICONOS_ROLLING_FRICTION_2D_ONECONTACT_ProjectionOnCone = 4001,
  /** Projection on cone with local iteration for 2D rolling friction */
  SICONOS_ROLLING_FRICTION_2D_ONECONTACT_ProjectionOnConeWithLocalIteration = 4002,
  /** Enumeration of the cases of one contact for 2D rolling friction (exact) */
  SICONOS_ROLLING_FRICTION_2D_ONECONTACT_ENUM = 4003,

  /* -----------------------------------------------------------------------
   * Global Rolling Friction solvers (IDs 5000-5999)
   * -----------------------------------------------------------------------
   * These solvers work on the global (aggregated) formulation of the
   * rolling friction contact problem.
   */
  /** Non-smooth Gauss-Seidel for 3D rolling friction, global formulation */
  SICONOS_GLOBAL_ROLLING_FRICTION_3D_NSGS_WR = 5000,
  /** Interior Point Method for 3D rolling friction */
  SICONOS_GLOBAL_ROLLING_FRICTION_3D_IPM = 5001
};

/* ===========================================================================
 * Solver Name Strings
 * ===========================================================================
 * Human-readable string representations of rolling friction solver IDs.
 */

#endif /* ROLLING_FRICTION_CONTACT_OPTIONS_H */
