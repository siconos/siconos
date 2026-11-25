
/* Siconos is a program dedicated to modeling, simulation and control
 * of non smooth dynamical systems.
 *
 * SPDX-FileCopyrightText: INRIA
 * SPDX-License-Identifier: Apache-2.0
 */

#include "GlobalFrictionContactProblem_as_ConvexQP.h"

#include "ConvexQP.h"                      // for ConvexQP
#include "GlobalFrictionContactProblem.h"  // for GlobalFrictionContactProblem
#include "SiconosBlas.h"                   // for cblas_dcopy
#include "projectionOnCone.h"              // for projectionOnCone
/* #define DEBUG_STDOUT */
/* #define DEBUG_MESSAGES */
#include "numerics_errors.h"
#include "siconos_debug.h"

void Projection_ConvexQP_GFC3D_DualCone(void* cqpIn, double* x, double* PX) {
  DEBUG_PRINT("Projection_ConvexQP_FC3D_Cylinder(void *cqpIn, double *x, double *PX)\n")

  ConvexQP* cqp = (ConvexQP*)cqpIn;
  GlobalFrictionContactProblem_as_ConvexQP* pb =
      (GlobalFrictionContactProblem_as_ConvexQP*)cqp->env;
  GlobalFrictionContactProblem* gfc3d = pb->gfc3d;

  // globalFrictionContactProblem_display(gfc3d);

  int contact = 0;
  int nLocal = gfc3d->dimension;
  int n = gfc3d->numberOfContacts * nLocal;

  cblas_dcopy(n, x, 1, PX, 1);
  for (contact = 0; contact < gfc3d->numberOfContacts; ++contact) {
    projectionOnCone(&PX[contact * nLocal], 1.0 / gfc3d->mu[contact]);
  }
}
