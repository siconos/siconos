/* Siconos is a program dedicated to modeling, simulation and control
 * of non smooth dynamical systems.
 *
 * SPDX-FileCopyrightText: INRIA
 * SPDX-License-Identifier: Apache-2.0
 */

/*!\file rolling_friction_3d_short_names.h
 * \brief Short macro aliases for long RollingFrictionContact enum names
 *
 * This header provides shorter names for frequently used enum values.
 * It is automatically included by RollingFrictionContact headers.
 *
 * Usage: Instead of SICONOS_ROLLING_FRICTION_3D_NSGS, use RFC3D_NSGS
 *        Instead of SICONOS_GLOBAL_ROLLING_FRICTION_3D_NSGS_WR, use GRFC3D_NSGS_WR
 *
 * The short names follow the pattern:
 *   RFC3D_*   = 3D Rolling Friction Contact solvers (local formulation)
 *   RFC2D_*   = 2D Rolling Friction Contact solvers
 *   GRFC3D_*  = Global Rolling Friction Contact 3D solvers
 *   GRFC2D_*  = Global Rolling Friction Contact 2D solvers
 */

#ifndef RFC3D_SHORT_NAMES_H
#define RFC3D_SHORT_NAMES_H

#include "RollingFrictionContact_options.h"

/* ===========================================================================
 * 3D Rolling Friction Contact Solvers - Local Formulation (RFC3D_*)
 * =========================================================================== */
#define RFC3D_NSGS SICONOS_ROLLING_FRICTION_3D_NSGS
#define RFC3D_OC_PROJ SICONOS_ROLLING_FRICTION_3D_ONECONTACT_ProjectionOnCone
#define RFC3D_OC_PROJ_LI \
  SICONOS_ROLLING_FRICTION_3D_ONECONTACT_ProjectionOnConeWithLocalIteration
#define RFC3D_ADMM SICONOS_ROLLING_FRICTION_3D_ADMM

/* ===========================================================================
 * 2D Rolling Friction Contact Solvers (RFC2D_*)
 * =========================================================================== */
#define RFC2D_NSGS SICONOS_ROLLING_FRICTION_2D_NSGS
#define RFC2D_OC_PROJ SICONOS_ROLLING_FRICTION_2D_ONECONTACT_ProjectionOnCone
#define RFC2D_OC_PROJ_LI \
  SICONOS_ROLLING_FRICTION_2D_ONECONTACT_ProjectionOnConeWithLocalIteration
#define RFC2D_OC_ENUM SICONOS_ROLLING_FRICTION_2D_ONECONTACT_ENUM

/* ===========================================================================
 * Global Rolling Friction Contact 3D Solvers (GRFC3D_*)
 * =========================================================================== */
#define GRFC3D_NSGS_WR SICONOS_GLOBAL_ROLLING_FRICTION_3D_NSGS_WR
#define GRFC3D_IPM SICONOS_GLOBAL_ROLLING_FRICTION_3D_IPM

#endif /* RFC3D_SHORT_NAMES_H */
