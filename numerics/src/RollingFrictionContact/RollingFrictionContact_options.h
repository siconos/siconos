/* Siconos is a program dedicated to modeling, simulation and control
 * of non smooth dynamical systems.
 *
 * Copyright 2024 INRIA.
 *
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 *
 * http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 */

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
  /** ADMM for 3D rolling friction, local formulation */
  SICONOS_ROLLING_FRICTION_3D_ADMM = 3003,

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

/* ===========================================================================
 * Short Name Aliases (Optional)
 * ===========================================================================
 * For convenience, shorter macro aliases are available by including
 * rfc3d_short_names.h. These provide shorter names like RFC3D_NSGS instead of
 * SICONOS_ROLLING_FRICTION_3D_NSGS.
 *
 * Usage:
 *   #include "rfc3d_short_names.h"
 *   // Or use the long names directly (backward compatible)
 */

#endif /* ROLLING_FRICTION_CONTACT_OPTIONS_H */
